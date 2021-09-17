source("../simulator/simulator.R")

#### tests #####

counts_maynard <- Matrix(as.matrix(fread("~/ma/data/Maynard/X_tpm.csv")), sparse = T)
genes <- fread("~/ma/data/Maynard/var.csv")
cells <- fread("~/ma/data/Maynard/obs.csv")
cellnames <- cells$Run
genenames <- genes$symbol
dimnames(counts_maynard)<-list(cellnames, genenames)
counts_maynard <- t(counts_maynard)

annotation_maynard <- fread("~/ma/data/Maynard/annotation_obs.csv")[,c("Run","cell_type")]
colnames(annotation_maynard) <- c("ID", "cell_type")

ds_m <- dataset(annotation_maynard, counts_maynard, name = "Maynard", count_type = "TPM", 
                whitelist = c("B cell","Macrophage","NK cell","T cell CD8","T cell CD4","T cell regulatory","Monocyte conventional","Monocyte non-conventional"))

ds_m2 <- dataset(annotation_maynard, counts_maynard, name = "Maynard", count_type = "TPM")

# annotation_travaglini <- fread("~/ma/data/Travaglini/obs_extended.csv")
# counts_travaglini <- "~/ma/data/Travaglini/Travaglini_Krasnow_2020_Lung_SS2.h5ad"
# 
# ercc_cols <- grep("ERCC-",colnames(annotation_travaglini))
# meta_red <- data.frame(annotation_travaglini)[ercc_cols]
# meta_red$ercc_sum <- apply(meta_red,1,sum)
# meta_red$cell_type <-annotation_travaglini$free_annotation
# meta_red$total_reads <- annotation_travaglini$nReads
# meta_red$genes <- annotation_travaglini$nGene
# meta_red$ID <- annotation_travaglini$index
# meta_red <- as.data.table(meta_red)
# 
# ds_t <- dataset_h5ad(meta_red, counts_travaglini, name = "Travaglini", spike_in_col = "ercc_sum")
# 
# db <- database(list(ds_m, ds_t))

db <- database(list(ds_m))

sim<-simulate_bulk(ds, 
                   scenario = "random", 
                   scaling_factor="NONE",
                   ncores=4, 
                   ncells=1000, 
                   nsamples=100) 

# sim$pseudo_bulk %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(-c(gene), names_to = "samples", values_to = "counts") %>%
#   ggplot(aes(x=samples,y=gene,fill=counts))+
#   geom_raster()+
#   scale_fill_viridis_c()