# Sekine-et-al.-2020

## 1. Download data
scRNA-seq data were obtained from the public 10× genomics application note: _A New Way of Exploring Immunity_(www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexedantigen-recognition-to-immune-repertoire-and-phenotype/). 

* [donor 1](https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor1)  
`curl -O https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor1/vdj_v1_hs_aggregated_donor1_filtered_feature_bc_matrix.tar.gz`
* [donor 4](https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.2/vdj_v1_hs_aggregated_donor4)  
`curl -O https://cf.10xgenomics.com/samples/cell-vdj/3.0.2/vdj_v1_hs_aggregated_donor4/vdj_v1_hs_aggregated_donor4_filtered_feature_bc_matrix.tar.gz`

## 2. Normalize and integrate
CD8+ T cell gene expression data from both healthy donors was then normalized with SCTransform and integrated with Seurat v3 using the script [integrate_d1_d4.R](./integrate_d1_d4.R)

## 3. DEG
Tox+ (1909) and Tcf+ (9430) cells were identified, and differentially expressed genes (DEGs) were determined using DESeq2 using the script [DEG_Tox_Tcf7.R](./DEG_Tox_Tcf7.R)
