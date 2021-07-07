maynard <- fread("../../Data/Maynard/read_count_per_cell_type.csv")
colnames(maynard) <- c("type","mRNA")
#maynard$mRNA <- range01(maynard$mRNA)+0.1
maynard$type <- c("Ciliated", "Goblet", "Epithelial cell other C4", "Mesothelial", "Alevolar cell type 2", "Plasma cells", "Epithelial cell other C5", "Club cells", "T cells dividing", "Fibroblast adventitial", "Epithelial cell other C3", "Basal", "Endothelial cells L", "Ionocyte", "conventional Dendritic cells 2", "Epithelial cell other C2", "Dendritic cells", "other C21", "Epithelial cell other C1", "Mast cells", "conventional Dendritic cells 1", "Endothelial cell", "Alevolar cell type 1", "Macrophages", "plasmacytoid Dendritic cells", "T regulatory cells", "Smooth muscle cells", "Pericyte", "Fibroblast alevolar", "B cells", "other C33", "other C29", "NK cells", "Monocyte non-conventional", "T cells CD8", "other C32", "T cells CD4", "Monocytes", "Dendritic cells Langerhans", "other C19", "other C17", "Epithelial cell dividing")

quantiseq <- fread("../../Data/quanitseq.tsv")
colnames(quantiseq) <- c("type","mRNA")
#quantiseq$mRNA <- range01(quantiseq$mRNA)+0.1
quantiseq$type <- c("B cells","Macrophages","MacrophagesM2","Monocytes","Neutrophils","NK cells","T cells CD4","T cells CD8","T regulatory cells","Dendritic cells")

load("../../Data/epic.rda")
epic <- data.frame(mRNA_cell_default)
epic$type <- c("B cells","Macrophages","Monocytes","Neutrophils","NK cells" ,"T cells","T cells CD4","T cells CD8","T helper cells","T regulatory cells","otherCells","default")
colnames(epic)<-c("mRNA","type")
#epic$mRNA <- range01(epic$mRNA)+0.1

miltenyi <- data.frame(type=c("B cells", "Hematopoietic progenitor cells", "Monocytes", "conventional Dendritic cells", "plasmacytoid Dendritic cells","NK cells", "NKT cells", "T cells", "T cells CD4", "T cells CD8", "T regulatory cells"), mRNA=c(0.4, 0.5, 1, 2.5, 0.7, 0.3, 0.6, 1, 0.6, 0.4, 0.3))
#miltenyi$mRNA <- range01(miltenyi$mRNA)+0.1

meta_runs <- fread("../../Data/Vento-Tormo/filereport_read_run_PRJEB25794_tsv.txt")[,c("library_name","read_count")]
meta_cells <- fread("../../Data/Vento-Tormo/meta_ss2.txt")[,c("V1","annotation")]
meta_cells$ID <- gsub("#", "_", meta_cells$V1)
meta_runs$library_name <- gsub("_p", "", meta_runs$library_name)
meta <- merge(meta_cells,meta_runs, by.x="ID", by.y = "library_name")
colnames(meta) <- c("id","id2","cell_type","n_reads")
meta[,cell_count:=.N, by="cell_type"]
vento <- meta
vento_tormo <- vento[, mean(n_reads), by=cell_type]
colnames(vento_tormo) <- c("type","mRNA")
vento_tormo$type <- c("villous cytotrophoblast", "decidual Macrophages 1", "decidual Macrophages 2", "T cells", "conventional Dendritic cells 1", "decidual NK cells 2", "decidual NK cells 1", "decidual NK cells 3" ,"NK cells proliferating", "stromal cells 2", "stromal cells 1", "stromal cells 3","Endothelial cells (m)", "epithelial glandular cells 2", "innate lymphocyte cells", "dP2", "NK cells CD56bright CD16+/-", "Plasma cells", "Macrophages", "Hoffbauer cells", "Granulocytes", "Fibroblasts 1", "extravillous trophoblast", "syncytiotrophoblast", "dP1", "Endothelial cells L", "epithelial glandular cells 1", "Fibroblasts 2", "Endothelial cells (f)", "NK cells CD16 -", "Monocytes", "conventional Dendritic cells 2")

cell_types <- c("Monocytes", "NK cells", "T cells CD8 central memory", "T cells CD4","T cells CD8","B Naive" ,"T cells CD4 central memory","MAIT" ,"T gd Vd2" ,"Neutrophils" ,"T gd non-Vd2" ,"Basophils LD", "Monocytes NC+I" ,"B cells"   ,"conventional Dendritic cells"   ,"plasmacytoid Dendritic cells" ,"Plasma cells") 
fact <- c(22824.322779,  21456.905980,  17827.199342,  14262.064545,  10660.952413,  24162.648934,  15682.087676,  13365.625992,  28128.898864,  9546.742675,  25008.840273,  8610.360645,  22344.663722,  20837.569855,  57322.183941,  28488.194770,  325800.989149)
monaco <- data.frame(type=cell_types, mRNA=fact)

df <- fread("../../Data/Hao (CITEseq-PBMC)/census_meta.tsv")
census_mean_counts <-df[,mean(census_counts),by=celltype.l2]
census_mean_counts$source <- "census (mean)"
colnames(census_mean_counts) <- c("type","mRNA","source")
census_mean_counts$type <- c("T cells CD8 effector memory","B cells","T cells CD4 central memory","Monocytes CD14", "T cells CD4", "Platelet", "T cells CD4 effector memory", "T cells CD8", "NK cells", "Monocytes CD16", "B cells intermediate", "B cells memory", "conventional Dendritic cells 2", "T cells CD8 central memory", "HSPC", "T cells gamma delta", "NK cells CD56bright CD16+/-", "T cells CD4 cytotoxic activity", "Doublet", "plasmacytoid Dendritic cells", "MAI T cells", "T regulatory cells", "Plasma cells", "T cells CD4 proliferating", "ASDC", "NK cells proliferating", "innate lymphocyte cells", "T cells double negative", "conventional Dendritic cells 1", "Eryth", "T cells CD8 proliferating")
