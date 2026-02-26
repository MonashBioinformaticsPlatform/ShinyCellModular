
# Library load

library(ShinyCell) #devtools::install_github("SGDDNB/ShinyCell")
library(Seurat)
library(Signac)
library(dplyr)

setwd("/Users/lper0012/Library/CloudStorage/Dropbox/BioPlatform/ShinyCellPlus_devel/ShinyCellPlus/RProjShinyCellPlus")

##########################################################################################################
# Seurat and check Key names
##########################################################################################################

cnts<-LoadSeuratRds("../../input_multiomics/seurat_UMAP_clustered_integration_peaklinks_motifs_regulons.rds")

DefaultAssay(cnts)<-'RNA'
Idents(cnts)<-cnts$integrated_snn_res.1

# See what assays and keys you have
Assays(cnts)

for(a in Assays(cnts)) {
  cat(a, "key:", Key(cnts[[a]]), "\n")
}

# add key names to the assays that dont have

for(a in Assays(cnts)) {
  if (is.null(Key(cnts[[a]])) || nchar(Key(cnts[[a]])) == 0) {
    Key(cnts[[a]]) <- paste0(tolower(a), "_")
  }
}

# Check again

for(a in Assays(cnts)) {
  cat(a, "key:", Key(cnts[[a]]), "\n")
}

# find most variable features

cnts<-FindVariableFeatures(cnts)


##########################################################################################################
# Create market genes table
##########################################################################################################

# # install.packages("devtools")
#devtools::install_github("immunogenomics/presto")
library(presto)
DefaultAssay(cnts)<-'RNA'
resolutions<-colnames(cnts@meta.data)[grep("res.",colnames(cnts@meta.data))]

markers_list<-data.frame()

for (res in resolutions){
  Idents(cnts)<- cnts@meta.data[[res]]
  counts_cnts <- GetAssayData(cnts, layer = "data")
  cluster_labels <- cnts@meta.data[[res]]
  markers<- presto::wilcoxauc(counts_cnts, cluster_labels) |> as.data.frame() |> mutate(annotation=res)
    markers_list<-rbind(markers_list, markers)
}

library(arrow)
write_parquet(markers_list,"Files_ShinyCell/markergenes_lists.parquet")


##########################################################################################################
# add a reduction for 3D plotting
##########################################################################################################

cnts <- RunUMAP(
  cnts,
  reduction = "pca",  # or "pca"
  dims = 1:30,
  n.components = 3,
  reduction.name = "pca_umap3d"
)





##########################################################################################################
# Create the configuration files for a ShinyCell app (For only RNA assays)
##########################################################################################################
library(Seurat)
library(ShinyCell)

scConf = createConfig(cnts)
makeShinyApp( cnts, scConf, gene.mapping = TRUE, 
              shiny.title = "ShinyCell Quick Start",
              shiny.dir = "Files_ShinyCell") 

##########################################################################################################
# Create market genes table (For Multiomics or Spatial)
##########################################################################################################

library(Seurat)
library(ShinyCell2)

scConf <- createConfig(cnts)
makeShinyFiles(cnts,scConf, shiny.prefix="sc1", shiny.dir="Files_ShinyCell2/")
makeShinyCodes(shiny.title = "ShinyCell2 Pando Processed multiomics", shiny.prefix="sc1",
               shiny.dir="Files_ShinyCell2/")



##########################################################################################################
# Function that creates the ShinyCellPlus App
##########################################################################################################

source("../functions//useShinyCellPlus.R")


useShinyCellPlus(
    shiny.dir="Files_ShinyCell/",
    shinycellplus.dir.src="~/Dropbox/BioPlatform/ShinyCellPlus_devel/ShinyCellPlus/",
    rsconnect.deploy = FALSE,
    data_type = "",
    enabled_tabs = c("cellinfo_cellinfo","cellinfo_geneexpr","cellinfo3D_cellinfo3D","cellinfo3D_geneexpr3D","genecoex","violin_boxplot","proportions","bubble_heatmap","pseudobulk","testing"),
    overwrite_code = TRUE,
    app_title='Testing'
)



##########################################################################################################
# Adding Pesudobulking:
##########################################################################################################

library(hdf5r)
library(Matrix)


counts <- GetAssayData(cnts, assay = "RNA", layer= "counts")

dirname<- "Files_ShinyCell/"
file.h5 <- H5File$new(paste0(dirname,"sc1counts.h5"), mode = "w")
grp <- file.h5$create_group("counts")

counts_dense <- as.matrix(counts)
storage.mode(counts_dense) <- "integer"

grp$create_dataset(
  "counts",
  robj       = counts_dense,
  dtype      = h5types$H5T_STD_I32LE,
  chunk_dims = c(1, ncol(counts_dense)),
  gzip_level = 4
)

file.h5$close_all()

# ##########################################################################################################
# Function that creates a template to create a new module script
##########################################################################################################
# this is buggy doesnt no create sc1 prefixes corrrectly
createSCPModuleTemplate(module_dir = "Files_ShinyCell/code/modules/",
                  tab_id = "pseudobulk",
                  module_name="pseudobulk",
                  tab_title = "Pseudobulking")


# ##########################################################################################################
# prepare files wrapper
##########################################################################################################

prepShinyCellPlus(seurat_rds = "../../input_multiomics/seurat_UMAP_clustered_integration_peaklinks_motifs_regulons.rds",out_dir = "testing_data", do_umap3d = TRUE,  do_markers= TRUE,   markers_res_pattern = "wsnn_res")


















































































































library(dplyr)    # alternatively, this also loads %>%
library(Seurat)
library(hdf5r)
library(fs)
library(scCustomize)
library(qs)
library(clustree)
library(shiny)
library(SeuratDisk)
library(ShinyCell)
library(rhdf5)
setwd("~/data/tasks/margo.montandon/Analysis/SCENIC_MARCH2024/forSeurat/")
cnts<-LoadSeuratRds("sobj.merged_withLinkPeaks.Rds")
DefaultAssay(cnts)<-'RNA'
Idents(cnts)<-cnts$wsnn_res.0.1
cnts<-FindVariableFeatures(cnts)
scConf = createConfig(cnts)
makeShinyApp( cnts, scConf, gene.mapping = TRUE, shiny.title = "ShinyCell Quick Start", shiny.dir="ShinyCell_Zebrafish/configFiles_RNA") 


DefaultAssay(cnts)<-'ATAC'
Idents(cnts)<-cnts$wsnn_res.0.1
cnts<-FindVariableFeatures(cnts)
scConf = createConfig(cnts)
makeShinyApp( cnts, scConf, gex.assay = 'ATAC', gene.mapping = TRUE, shiny.title = "ShinyCell Quick Start", shiny.dir="ShinyCell_Zebrafish/configFiles_ATAC") 
linked_genes_ATAC <- Links(cnts)

saveRDS(linked_genes_ATAC,"ShinyCell_Zebrafish/configFiles_ATAC/sc1Links.rds")


##########################################################################################################################################
#from seurat object
unique(Links(cnts@assays$ATAC)$gene)



setwd("/mnt/ceph/mbp/servers/bioinformatics-platform/home/lper0012/tasks/margo.montandon/Analysis/SCENIC_MARCH2024/forSeurat/")


# RNA
dir_inputs="ShinyCell_Zebrafish/configFiles_RNA/"

sc1conf_RNA = readRDS(paste0(dir_inputs,"sc1conf.rds"))
sc1def_RNA  = readRDS(paste0(dir_inputs,"sc1def.rds"))
sc1gene_RNA = readRDS(paste0(dir_inputs,"sc1gene.rds"))
sc1meta_RNA = readRDS(paste0(dir_inputs,"sc1meta.rds"))

h5file_RNA <- H5File$new("ShinyCell_Zebrafish/configFiles_RNA/sc1gexpr.h5", mode = "r") 
h5data_RNA <- h5file_RNA[["grp"]][["data"]] 


# ATAC
dir_inputs="ShinyCell_Zebrafish/configFiles_ATAC/"

sc1conf_ATAC = readRDS(paste0(dir_inputs,"sc1conf.rds"))
sc1def_ATAC  = readRDS(paste0(dir_inputs,"sc1def.rds"))
sc1gene_ATAC= readRDS(paste0(dir_inputs,"sc1gene.rds"))
sc1meta_ATAC = readRDS(paste0(dir_inputs,"sc1meta.rds"))
sc1link_ATAC = readRDS(paste0(dir_inputs,"sc1Links.rds"))

h5file_ATAC <- H5File$new("ShinyCell_Zebrafish/configFiles_ATAC/sc1gexpr.h5", mode = "r") 
h5data_ATAC <- h5file_ATAC[["grp"]][["data"]] 

####################################################################################################






##############################################################
# what parts of the object I need

library(Signac)
library(Seurat)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)


# find the objects that have the links
sobj.merged<-LoadSeuratRds("~/data/tasks/margo.montandon/Analysis/SCENIC_MARCH2024/forSeurat/sobj.merged_withLinkPeaks.Rds")

# linkplot

# Gviz

CoveragePlot(sobj.merged, 
             region = c("19-22216256-22217163"),
             features = "nfatc1",
             expression.assay = "RNA" , 
             extend.upstream = 1000000,
             extend.downstream = 100000,
             links = TRUE,
             split.by = "sample")




h5data_ATAC 
> loadShinyCellATAC <- function(h5_path) {
  +     h5file <- H5File$new(h5_path, mode = "r")
  +     h5grp <- h5file[["grp"]]
  +     list(
    +         rownames = h5grp[["rownames"]][],
    +         colnames = h5grp[["colnames"]][],
    +         data = h5grp[["data"]]
    +     )
  + }


library(GenomicRanges)
library(Matrix)
library(ggplot2)
library(dplyr)

# sc1link_ATAC: GRanges object with gene name metadata
# atac_counts: sparse matrix of peak accessibility (peaks x cells)
# sc1meta: data.frame with cell annotations
# group_by: metadata column to summarize by

library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(Matrix)
library(hdf5r)

plotGeneAccessibility <- function(gene,
                                  sc1link_ATAC,
                                  h5data_ATAC,
                                  sc1meta,
                                  group_by = "cluster") {
  # Step 1: Get peaks linked to gene
  linked_peaks <- sc1link_ATAC[which(mcols(sc1link_ATAC)$gene == gene)]
  
  if (length(linked_peaks) == 0) {
    warning(paste("No peaks linked to gene:", gene))
    return(NULL)
  }
  
  # Step 2: Get peak names that match rownames in h5 data
  peak_names <- names(linked_peaks)
  h5_peak_names <- h5data_ATAC[["rownames"]]$read()
  peak_idx <- which(h5_peak_names %in% peak_names)
  
  if (length(peak_idx) == 0) {
    warning(paste("No matching peaks found in sc1gexpr.h5 for gene:", gene))
    return(NULL)
  }
  
  # Step 3: Read only necessary accessibility data (peaks x all cells)
  counts_mat <- h5data_ATAC[["data"]][peak_idx, ]
  
  # Get column (cell) names
  cell_names <- h5data_ATAC[["colnames"]]$read()
  
  # Match metadata
  meta <- sc1meta[match(cell_names, sc1meta$cell), ]
  
  if (!group_by %in% colnames(meta)) {
    stop(paste("Grouping column", group_by, "not found in metadata"))
  }
  
  # Assign group labels
  meta$group <- meta[[group_by]]
  valid_cells <- !is.na(meta$group)
  meta <- meta[valid_cells, ]
  counts_mat <- counts_mat[, valid_cells, drop = FALSE]
  
  # Compute group-wise mean accessibility across peaks
  df_access <- data.frame(cell = meta$cell, group = meta$group)
  df_access$access <- colMeans(counts_mat)
  
  df_plot <- df_access %>%
    group_by(group) %>%
    summarise(accessibility = mean(access, na.rm = TRUE)) %>%
    ungroup()
  
  # Plot
  ggplot(df_plot, aes(x = group, y = accessibility, fill = group)) +
    geom_col() +
    theme_minimal() +
    labs(title = paste("Accessibility for", gene),
         x = group_by,
         y = "Average accessibility (linked peaks)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


plotGeneAccessibility(
  gene = "nfatc1",
  sc1link_ATAC = sc1link_ATAC,
  h5data_ATAC = h5data_ATAC,
  sc1meta = sc1meta,
  group_by = "celltype"
)
