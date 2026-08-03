# ShinyCellModular → ShinyProxy portable deployment — design notes

Status as of 2026-07-17: **discussion only, no code changed yet.** This file is a
checkpoint to resume the conversation. Everything below is decided/agreed unless
marked "OPEN".

## Goal

Laura wants to deploy ShinyCellModular apps behind ShinyProxy. Tyron will run his
own ShinyProxy instance on a Nectar VM. Requirement given by Tyron/infra (not
Laura's own design): capture a `projectid` from the URL query string
(`{url}?projectid=bla-foo-hex-code`) and resolve it against a shared, read-only
mounted directory of projects, per this Posit article:
https://shiny.posit.co/r/articles/build/client-data/

## What already existed before this conversation

- `Dockerfile` (repo root): builds a single shared image (`lperlaza/shinycellmodular`
  on Docker Hub), installs all deps + the package itself. Current `CMD` runs
  whatever app is at `Sys.getenv('SCMODULAR_DATA_DIR', '/data')` — i.e. today's
  model is "one bundle mounted per container," not the new shared-mount model.
- `.github/workflows/docker-image.yaml`: builds/pushes that image on changes to
  `Dockerfile` or itself.
- `R/useShinyCellModular.R`: generates a **self-contained, local-dev** `app.R` into
  `out_dir` from an internal string template (`template_app`), baking
  `app_title`/`navbar_css`/`shinyCellModularVersion`/`dir_inputs`/`assays`/
  `enabled_tabs` in as **literal values** via `gsub("__X__", ...)`. No config file
  needed to run it locally. Also copies selected `modules/*.R` into `out_dir`.
- `R/makeShinyCellModularPortable.R`: takes an existing `out_dir` (built by
  `useShinyCellModular()`) and converts it into a deployable `bundle_dir`:
  copies `modules/` and data folders, writes `app_config.yml` (title, assays,
  enabled_tabs, navbar_css, scm_version), and rewrites `out_dir`'s app.R into a
  skeleton that reads those same values from `app_config.yml` at runtime instead
  of from baked literals. Values are **extracted** from `out_dir/app.R` via
  `source(app_r_path, local = source_env)` (deliberately not regex — this part
  is fine and was already correct). The **rewrite** step uses `sub()` regex
  matching the literal shape of the baked assignment block — this part *is*
  fragile (breaks if that block's text shape changes) and is the one thing we
  discussed fixing (see "Continuity fix" below).
- `R/templateShinyCellModular.R`: scaffolds a new module/tab file.
- `R/prepShinyCellModular.R`: Seurat object → ShinyCell config + data files
  (RNA/ATAC), unrelated to this deployment discussion.
- This is a single-dataset-per-container model end to end: one ShinyProxy
  container == one bundle_dir == one dataset.

## Architecture agreed for the ShinyProxy/Nectar/Ceph deployment

**Rejected**: one ShinyProxy `spec:` entry per dataset (my first sketch). Superseded
once Tyron's actual spec (query-string + shared mount) came up — that's strictly
better for a growing, Ceph-backed, multi-dataset setup.

**Agreed shape**:

```
Nectar instance running ShinyProxy
├── application.yml            — ONE spec (not one per dataset)
│     container-image: lperlaza/shinycellmodular:latest   (same image, unchanged)
│     container-volumes: ["<ceph-path>:/projects:ro"]     (whole shared tree, read-only)
│     container-env: PROJECTS_BASE=/projects
│
└── Ceph-backed shared tree (mounted at /projects inside every container)
    ├── _app/
    │   └── app.R              — ONE shared app.R, MOUNTED not baked into image
    │                             (reserved folder name `_app`, can't collide with
    │                             a real projectid since projectid is validated
    │                             against ^[A-Za-z0-9_-]+$ and folder must contain
    │                             app_config.yml to be treated as a project)
    ├── bla-foo-hex-code/       — one folder per project, projectid = folder name
    │   ├── app_config.yml
    │   ├── modules/            — CAN be custom per project (see decision below)
    │   └── RNA/ or ATAC/...    — data
    └── another-project-id/
        └── ...
```

Users reach a project via `https://<shinyproxy-host>/app_direct/<spec-id>/?projectid=bla-foo-hex-code`
(note: needs the slash before `?`, and ShinyProxy >= 3.0.2 for query-string
passthrough on `app_direct`/`app` — not yet confirmed which version Tyron runs,
should check before relying on this).

### Decisions locked in along the way

1. **Modules can be custom per project** (not a fixed set baked into the image).
   Consequence: module sourcing / `register_tab()` can't happen once at
   container-boot time — it has to happen once real request/session info
   (which project) is available.
2. **One container per browser session** (ShinyProxy default), not
   pooled/shared. Simplifies memory reasoning (no multi-tenant cache needed) but
   does **not** remove the need to move data/module loading out of top-level
   app.R — the container's R process starts and sources app.R fully *before*
   any browser ever connects, so `projectid` genuinely isn't known yet at that
   point regardless of pooling.
3. **`projectid`-in-URL is being used as the de facto access control** ("unique
   url = unique identifier"). Flagged caveat, not re-litigated: this only works
   as real access control if folder names are unguessable and ShinyProxy's app
   list/landing page doesn't enumerate them; not a substitute for real auth if
   that matters for the data sensitivity involved.
4. Checked ShinyProxy's own "App Parameters" feature (SpEL-templated
   `container-volumes`/`container-env`, resolved *before* the container starts,
   which would have let app.R stay completely unchanged). **Ruled out**: its
   docs are explicit that "all possible values [must be] hard coded in the
   configuration file... it's not possible to have the user enter a text
   message" — deliberate, to prevent free user text flowing into a container
   mount path. Doesn't fit an arbitrary/growing set of project IDs without
   editing `application.yml` and restarting ShinyProxy per new project, which
   defeats the point of the Ceph/shared-mount design.
5. Landing page (type a projectid, press launch) **does not need to be
   R/Shiny/ShinyProxy-aware at all** — it can be a trivial static HTML page
   anywhere that just redirects to
   `.../app_direct/<spec-id>/?projectid=<encodeURIComponent(typed value)>`.
   The one place the typed string actually gets resolved/validated into a real
   path has to be app.R (or a custom broker service reinventing part of what
   ShinyProxy already does — considered and rejected as unnecessary extra risk
   for no real gain).
6. **app.R must NOT be baked into the Docker image.** Laura's explicit
   requirement (matches her very first message: "mount ... even the app
   skeleton"): as she adds SPATIAL/CROPseq/MULTI support she'll be editing
   app.R's logic frequently and does not want an image rebuild+push+redeploy
   cycle just to ship an app.R change. Resolution: mount it too, at a fixed
   reserved path (`PROJECTS_BASE/_app/app.R`) alongside the projects tree;
   image `CMD` always points at that fixed path. Updating app.R = edit the
   file on shared storage; next session (fresh container) picks it up
   immediately, no rebuild needed. Image rebuild only needed when a *new
   system library / R package* is required (e.g. a Spatial package with new
   deps), not for app-logic changes.

## What the shared ShinyProxy app.R needs to contain (proposed shape)

Confirmed against Shiny's own docs (`session$clientData$url_search` is
**server-side only**, not available in `ui`; `parseQueryString()` is the
parsing helper) and against community reports that `ui <- function(request)`
reading `request$QUERY_STRING` is the standard way to get the *same* query
string before any session exists (needed because `ui` can't see
`session$clientData`).

Sections:

1. **Libraries + static globals** — unchanged from today's template: all the
   `library(...)` calls, `cList`, `pList`/`pList2`/`pList3`, `sList`, `lList`,
   `g_legend()`, `sctheme()`. None of these depend on `projectid`, stay at true
   top level.
2. **`resolve_project(query_string)`** helper — parses `projectid` out of
   either `request$QUERY_STRING` (called from `ui()`) or
   `session$clientData$url_search` (called from `server()`), validates it
   against `^[A-Za-z0-9_-]+$` (reject anything else — this is user input
   feeding a filesystem path), builds
   `project_dir <- file.path(Sys.getenv("PROJECTS_BASE","/projects"), pid)`,
   confirms `app_config.yml` exists there, reads it via `yaml::read_yaml()`.
3. **`load_project_data(project_dir, app_config)`** — today's top-level
   `if (file.exists(file.path(dir_inputs,"RNA")))` / `"ATAC"` blocks, moved
   into a function returning a named list (`sc1conf`, `sc1meta`, `sc1gene`,
   `sc1def`, the `_atac` equivalents, `markers_list`, etc.). This is the one
   place future SPATIAL/CROPseq/MULTI branches get added.
4. **`build_tab_registry(modules_dir)`** — replaces the current global
   `tab_registry <<-` pattern. Builds a fresh **local** list/env each call,
   defines `register_tab()` as a closure over it, sources every `.R` under
   `modules_dir`. Needs to be local (not global `<<-`) because different
   projects can have different custom modules, and a browser reconnect can spin
   up a second Shiny session inside the same still-running container.
5. **`ui <- function(request) { ... }`** — `resolve_project(request$QUERY_STRING)`,
   then `load_project_data()` (UI construction needs `sc1conf`/`sc1def` to
   populate widgets, per today's module contract) and `build_tab_registry()`,
   then builds the same `tabPanel`/`navbarPage`/footer structure as today,
   wrapped in `tryCatch` so an invalid/unknown `projectid` renders a plain
   "project not found" page instead of an R error dump.
6. **`server <- function(input, output, session) { ... }`** — same
   `resolve_project()` (independently re-parses the same query string from
   `session$clientData$url_search`), `load_project_data()`,
   `build_tab_registry()`, then the existing `args_to_pass`/`formals()`-matching
   dispatch loop, essentially unchanged in spirit from today.
7. **`shinyApp(ui, server)`** — same call, `ui` is just a function now.

**Dropped** from today's template for this mode: `get_app_dir()`/`modules_dir`
path-sniffing (modules_dir is always `file.path(project_dir, "modules")` now),
the global `tab_registry <<-` pattern.

**Accepted cost**: light metadata (`app_config.yml`, `sc1conf`/`sc1def`) gets
read twice per session — once in `ui()`, once in `server()` — because they're
different invocation contexts with no shared closure by default. Cheap (small
RDS/YAML reads, not the big h5/fragment data), so not worth the complexity of a
cache to avoid it. Could revisit with a projectid-keyed cache later if it ever
matters, since one-container-per-session means the cache would only ever hold
one entry per container's lifetime anyway.

**Not yet reconciled (OPEN, see below)**: how this `ui(request)`/`server(session)`
restructuring for the ShinyProxy multi-project app.R relates to `useShinyCellModular()`'s
template, given the very next decision below says `useShinyCellModular()`'s
generated app.R must stay a plain, self-contained, single-dataset, baked-values
script. These are two different shapes of app.R and the conversation had not
yet landed on how (or whether) they share code before compression happened.

## Continuity between `useShinyCellModular()` and `makeShinyCellModularPortable()`

Laura's explicit requirement, stated directly, superseding an earlier wrong
turn in this conversation:

> I want `useShinyCellModular` to be self contained, and run fine locally, not
> configuration file required — just a normal app.R. Once I have the shape of
> app.R and modules I want to make that portable.

**Rejected approach** (mine, wrong): making `useShinyCellModular()` itself
generate a config-driven (`app_config.yml`-reading) app.R from the start, to
eliminate the need for `makeShinyCellModularPortable()` to transform anything.
Explicitly vetoed — `useShinyCellModular()`'s output must stay self-contained
with baked literal values and zero config-file dependency, and Laura will keep
maintaining/extending that template directly (adding SPATIAL/CROPseq etc.
there) — that workflow is intentional, not a bug to fix.

**Agreed fix instead** — narrow, surgical, preserves everything above:

- Extraction of current values from `out_dir/app.R` stays exactly as it is
  today: `source(app_r_path, local = source_env)`. Already correct, already
  not regex, no change needed.
- The fragile part is the *rewrite* step — `sub()` matching the literal text
  shape of the baked assignment block
  (`'app_title <- "[^"]*"\n\nnavbar_css <- "[^"]*"\n\n...'`). This breaks if
  that block's shape changes at all.
- **Fix**: wrap that block in the template with two fixed sentinel comments
  that never need to change:
  ```r
  # ---- SCM_INSTANCE_VALUES_START ----
  app_title <- "__APP_TITLE__"
  navbar_css <- "__NAVBAR_CSS__"
  shinyCellModularVersion <- "__SCM_VERSION__"
  dir_inputs <- "__DIR_INPUTS__/"
  assays <- __ASSAYS__
  enabled_tabs <- __ENABLED_TABS__
  # ---- SCM_INSTANCE_VALUES_END ----
  ```
  `makeShinyCellModularPortable()` searches for those two literal marker
  strings (trivial, exact match — not fragile, since we control their exact
  text and they never change) and replaces everything between them with the
  config-driven equivalent block, regardless of what's inside or how it's
  formatted.
- **Ongoing maintenance contract**: whenever Laura adds a new baked
  instance-value to the template (e.g. a SPATIAL-specific path), she adds it
  inside the marker block in `useShinyCellModular()`, and adds the one matching
  line to `makeShinyCellModularPortable()`'s replacement block (pulling the
  equivalent value from `app_config.yml`). Small, mechanical, localized
  two-line edit each time — not a regex audit.
- This preserves the documented behavior that hand-edited `app.R`/`modules/`
  in `out_dir` get carried into the bundle unchanged (still transforming
  whatever text actually sits in `out_dir/app.R`, just anchored on markers
  instead of guessing at value patterns).

**Not yet done**: actually implementing this marker swap in
`R/useShinyCellModular.R` / `R/makeShinyCellModularPortable.R`. Last message
before compression asked "want me to implement this now?" — unanswered.

## Required steps for the ShinyProxy solution overall (revised recap)

1. New shared app.R for ShinyProxy mode — per "what needs to go in the app.R
   skeleton" section above. **Open question**: how much of this shares code
   with `useShinyCellModular()`'s template vs. being a genuinely separate
   third template. Not resolved yet.
2. ~~Bake app.R into the image~~ — **superseded**: mount it instead, at
   `PROJECTS_BASE/_app/app.R`, image `CMD` points at that fixed path.
3. A bundle builder for this mode: writes `<projectid>/app_config.yml`,
   `<projectid>/modules/`, `<projectid>/<data>` under the Ceph tree — no
   `app.R` per project needed (shared, mounted separately per #2). Could reuse/
   extend `makeShinyCellModularPortable()` or be a new lean sibling function —
   not decided.
4. ShinyProxy `application.yml`: one spec, `container-volumes` pointing at the
   Ceph path mounted to `/projects:ro`, `container-env: PROJECTS_BASE=/projects`.
   Confirm one-container-per-session (default, no pooling) — not yet written.
5. Validation in app.R: regex-validate `projectid`, confirm `app_config.yml`
   exists before doing anything else, clear error otherwise. Folded into
   `resolve_project()` above.
6. Landing page: static HTML, no R, text box building
   `.../app_direct/<spec-id>/?projectid=<typed>` — not yet built.

## Nothing has been implemented yet

No files have been changed in this repo as part of this discussion. Everything
above is design only. Suggested resumption order when picking this back up:

1. Decide/resolve the open question: does the ShinyProxy multi-project app.R
   share a template mechanism with `useShinyCellModular()`'s, or is it
   deliberately a separate, third thing? (This wasn't settled before
   compression — worth re-confirming with Laura first thing.)
2. Implement the marker-based continuity fix in
   `R/useShinyCellModular.R` / `R/makeShinyCellModularPortable.R` (smaller,
   self-contained, doesn't depend on resolving #1).
3. Build the ShinyProxy-mode bundle builder + shared app.R.
4. Update `Dockerfile`'s `CMD` to point at the mounted `_app/app.R` path
   (behind an env var or a second image variant — not yet decided which).
5. Write the `application.yml` template.
6. Write the static landing page.
