library(Seurat)
library(ggplot2)
library(patchwork)
library(glue)

options(future.globals.maxSize = 20000 * 1024^2)

read_raw_matrix <- function(donors, rds_fn) {
  print("checking matrix import")
  
  all_libraries <- list()
  
  if (!file.exists(rds_fn)) {
    all_libraries <- lapply(donors, function(donor) {
      print(glue("working on {donor}"))
      
      data_dir <- file.path("./data", donor)
      
      new_data <- Read10X(data.dir = file.path(data_dir, "filtered_feature_bc_matrix"))
      d_data <- CreateSeuratObject(project = donor, counts = new_data$`Gene Expression`)
      
      d_data[['Protein']] <- CreateAssayObject(counts = new_data$`Antibody Capture`)
      d_data[["donor.ID"]] <- factor(donor)
      
      # library IDs
      library_ID = read.csv(file.path(data_dir, "LibraryID.csv"), stringsAsFactors = FALSE)
      row.names(library_ID) <- library_ID$Barcode
      library_ID$LibraryID <- as.factor(library_ID$LibraryID)
      
      d_data <- AddMetaData(object = d_data, metadata = library_ID$LibraryID, col.name = 'library.ID')
      
      # user specific settings! copy/pasta from original code.
      
      ## mitochondrial fraction
      d_data[["percent.mt"]] <- PercentageFeatureSet(d_data, pattern = "^MT-")
      ## ribosomal fraction
      d_data[["percent.ribo"]] <- PercentageFeatureSet(d_data, pattern = "^RP[SL]")
      d_data <- subset(d_data,
                       subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA > 2600 & percent.mt < 7.5 & percent.ribo > 15)
      
      # Split libraries and Normalise with SCT transform
      lib.list <- SplitObject(d_data, split.by = "library.ID")
      lib.list <- lapply(lib.list, function(x) {
        x <- SCTransform(x)
        return(x)
      })
      return(lib.list)
    })
    
    all_libraries <- unlist(all_libraries)
    print(all_libraries)
    print("saving RDS file")
    saveRDS(all_libraries, file = rds_fn)
  } else {
    print("found RDS file")
    print("loading RDS file")
    all_libraries <- readRDS(rds_fn)
  }
  
  return(all_libraries)
}

integrate_datasets <- function(all_libraries, rds_fn) {
  print("these are the subsets we are working with:")
  chx <- c("vdj_v1_hs_unsorted_donor1_5gex_protein_1",
          "vdj_v1_hs_unsorted_donor1_5gex_protein_2",
          "vdj_v1_hs_sorted_donor1_5gex_protein_36",
          "vdj_v1_hs_sorted_donor1_5gex_protein_23",
          "vdj_v1_hs_sorted_donor1_5gex_protein_8",
          "vdj_v1_hs_sorted_donor4_5gex_protein_7",
          "vdj_v1_hs_sorted_donor4_5gex_protein_6",
          "vdj_v1_hs_sorted_donor4_5gex_protein_4",
          "vdj_v1_hs_unsorted_donor4_5gex_protein_1")

  all_libraries <- all_libraries[which(names(all_libraries) %in% chx)]

  print(names(all_libraries))

  print("starting integration")
  if (!file.exists(rds_fn)) {
    all_libraries_features <- SelectIntegrationFeatures(object.list = all_libraries,
                                                     nfeatures = 2000)
    
    all_libraries <- PrepSCTIntegration(object.list = all_libraries,
                                     anchor.features = all_libraries_features,
                                     verbose = TRUE)
    
    all_libraries <- lapply(X = all_libraries, FUN = RunPCA, verbose = FALSE, features = all_libraries_features)

    all_libraries_anchors <- FindIntegrationAnchors(object.list = all_libraries,
                                                 normalization.method = "SCT",
                                                 reduction = "rpca",
                                                 anchor.features = all_libraries_features,
                                                 verbose = TRUE)

    all_libraries_integ <- IntegrateData(anchorset = all_libraries_anchors,
                                      normalization.method = "SCT",
                                      verbose = TRUE)
    
    print("saving RDS file")
    saveRDS(all_libraries_integ, rds_fn)
  } else {
    print("found RDS file")
    print("loading RDS file")
    all_libraries_integ <- readRDS(rds_fn)
  }
  
  return(all_libraries_integ)
}

run_dim_reduction <- function(all_libraries_integ, rds_fn) {
  #all_libraries_integ <- FindVariableFeatures(all_libraries_integ)  

  print("starting PCA/UMAP")
  all_libraries_integ <- RunPCA(all_libraries_integ, verbose = TRUE)
  
  pca_elbow <- ElbowPlot(all_libraries_integ)
  ggsave("img/elbowplot.png", plot = pca_elbow, device = "png")
  
  all_libraries_integ <- RunUMAP(all_libraries_integ, dims = 1:30)
  
  saveRDS(all_libraries_integ, rds_fn)
}

# main code...
rds_dir <- "./rds"
if (!dir.exists(rds_dir)) {
  dir.create(rds_dir)
}

donors <- c(1, 4)
donors <- paste0("donor", donors)

all_libraries_fn <- glue("{rds_dir}/all_libraries_sctransformed.rds")
all_libraries_integ_fn <- glue("{rds_dir}/all_integ_sct_subset.rds")
umap_fn <- glue("{rds_dir}/all_integ_sct_subset_umap.rds")

if (!file.exists(umap_fn)) {
  all_libraries <- read_raw_matrix(donors, all_libraries_fn)
  all_libraries_integ <- integrate_datasets(all_libraries, all_libraries_integ_fn)

  rm(all_libraries)
  all_libraries_integ <- run_dim_reduction(all_libraries_integ, umap_fn)
} else {
  print("found umap rds...loading now")
  all_libraries_integ <- readRDS(umap_fn)
}


DefaultAssay(all_libraries_integ) <- "RNA"
all_libraries_integ <- NormalizeData(all_libraries_integ, verbose = TRUE)
# features <- c("TOX", "TCF7", "CCR7", "SELL", "IL7R", "S100A4", "GZMB", "GZMK", "GZMH")
features <- c("TOX", "TCF7", "CCR7", "IL7R", "SELL", "NELL2", "LEF1", "BACH2",
"MYC", "ID3", "PDCD1", "TIGIT", "CD244", "LAG3",
"GZMB", "GZMA", "PRF1", "CX3CR1", "TBX21", "EOMES",
"ZEB2", "PRDM1")

fp <- FeaturePlot(all_libraries_integ, features = features, pt.size = 0.2, ncol = 5)

ggsave(filename = "feature_sct_integ.svg", device = "svg", width = 12, height = 9)


