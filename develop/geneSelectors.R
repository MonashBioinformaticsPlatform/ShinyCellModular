geneSelectorsServer <- function(id, sc1conf, sc1gene, sc1def, fID) {
  moduleServer(id, function(input, output, session) {

    observe_helpers()

    optCrt <- "{ 
      option_create: function(data, escape) {
        return('<div class=\"create\"><strong>' + '</strong></div>');
      } 
    }"

    updateSelectizeInput(
      session, "sc1a1inp2",
      choices  = names(sc1gene),
      server   = TRUE,
      selected = sc1def$gene1,
      options  = list(
        maxOptions = 7,
        create     = TRUE,
        persist    = TRUE,
        render     = I(optCrt)
      )
    )

    updateSelectizeInput(
      session, "sc1a3inp1",
      choices  = names(sc1gene),
      server   = TRUE,
      selected = sc1def$gene1,
      options  = list(
        maxOptions = 7,
        create     = TRUE,
        persist    = TRUE,
        render     = I(optCrt)
      )
    )

    updateSelectizeInput(
      session, "sc1a3inp2",
      choices  = names(sc1gene),
      server   = TRUE,
      selected = sc1def$gene2,
      options  = list(
        maxOptions = 7,
        create     = TRUE,
        persist    = TRUE,
        render     = I(optCrt)
      )
    )

    updateSelectizeInput(
      session, "sc1b2inp1",
      choices  = names(sc1gene),
      server   = TRUE,
      selected = sc1def$gene1,
      options  = list(
        maxOptions = 7,
        create     = TRUE,
        persist    = TRUE,
        render     = I(optCrt)
      )
    )

    updateSelectizeInput(
      session, "sc1b2inp2",
      choices  = names(sc1gene),
      server   = TRUE,
      selected = sc1def$gene2,
      options  = list(
        maxOptions = 7,
        create     = TRUE,
        persist    = TRUE,
        render     = I(optCrt)
      )
    )

    updateSelectizeInput(
      session, "sc1c1inp2",
      server   = TRUE,
      choices  = c(sc1conf[is.na(fID)]$UI, names(sc1gene)),
      selected = sc1conf[is.na(fID)]$UI[1],
      options  = list(
        maxOptions = length(sc1conf[is.na(fID)]$UI) + 3,
        create     = TRUE,
        persist    = TRUE,
        render     = I(optCrt)
      )
    )
  })
}
