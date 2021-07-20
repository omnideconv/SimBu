source("/home/simulator/scripts/census.R")

load_maynard <- function(){
  maynard_reads_fractions <- fread("/home/Data/Maynard/read_count_per_cell_type.csv")
  colnames(maynard_reads_fractions) <- c("type","mRNA")
  #maynard_reads$type <- c("Ciliated", "Goblet", "Epithelial cell other C4", "Mesothelial", "Alevolar cell type 2", "Plasma cells", "Epithelial cell other C5", "Club cells", "T cells dividing", "Fibroblast adventitial", "Epithelial cell other C3", "Basal", "Endothelial cells L", "Ionocyte", "conventional Dendritic cells 2", "Epithelial cell other C2", "Dendritic cells", "other C21", "Epithelial cell other C1", "Mast cells", "conventional Dendritic cells 1", "Endothelial cell", "Alevolar cell type 1", "Macrophages", "plasmacytoid Dendritic cells", "T regulatory cells", "Smooth muscle cells", "Pericyte", "Fibroblast alevolar", "B cells", "other C33", "other C29", "NK cells", "Monocyte non-conventional", "T cells CD8", "other C32", "T cells CD4", "Monocytes", "Dendritic cells Langerhans", "other C19", "other C17", "Epithelial cell dividing")
  
  tpm <- fread("/home/Data/Maynard/datasets_annotated/maynard_2020_annotated_fine/X_tpm.csv")
  mat <- Matrix(as.matrix(tpm), sparse = T)
  rm(tpm)
  genes <- fread("/home/Data/Maynard/datasets_annotated/maynard_2020_annotated_fine/var.csv")
  cells <- fread("/home/Data/Maynard/datasets_annotated/maynard_2020_annotated_fine/obs.csv")
  cellnames <- cells$Run
  genenames <- genes$symbol
  dimnames(mat)<-list(cellnames, genenames)
  
  maynard_census <- census(t(mat), expr_threshold = 0.1,ncores=1, method="monocle")
  maynard_census <- data.frame(census=maynard_census, Run=names(maynard_census))
  
  meta_runs <- fread("/home/Data/Maynard/SRA_meta.txt")[,c("Run","AvgSpotLen","Bases")]
  meta_runs$n_reads <- meta_runs$Bases/meta_runs$AvgSpotLen
  meta_cells <- fread("/home/Data/Maynard/annotation_obs.csv")[,c("Run","cell_type")]
  maynard_reads <- merge(meta_cells, meta_runs, by.x = "Run", by.y="Run")
  
  maynard_both <- data.table(merge(maynard_census, maynard_reads, by.x="Run", by.y="Run"))
  maynard_census_fractions <- maynard_both[,mean(census), by=cell_type]
  colnames(maynard_census_fractions) <- c("type","frac")

  final <- merge(maynard_reads_fractions, maynard_census_fractions, by="type", sort=F)
  final$type <- c("Ciliated", "Goblet", "Epithelial cell other C4", "Mesothelial", "Alevolar cell type 2", "Plasma cells", "Epithelial cell other C5", "Club cells", "T cells dividing", "Fibroblast adventitial", "Epithelial cell other C3", "Basal", "Endothelial cells L", "Ionocyte", "conventional Dendritic cells 2", "Epithelial cell other C2", "Dendritic cells", "other C21", "Epithelial cell other C1", "Mast cells", "conventional Dendritic cells 1", "Endothelial cell", "Alevolar cell type 1", "Macrophages", "plasmacytoid Dendritic cells", "T regulatory cells", "Smooth muscle cells", "Pericyte", "Fibroblast alevolar", "B cells", "other C33", "other C29", "NK cells", "Monocyte non-conventional", "T cells CD8", "other C32", "T cells CD4", "Monocytes", "Dendritic cells Langerhans", "other C19", "other C17", "Epithelial cell dividing")
  
  maynard_census_fractions <- final[,c("type","frac")]
  maynard_census_fractions[maynard_census_fractions$type=="conventional Dendritic cells",]$frac<-mean(maynard_census_fractions[maynard_census_fractions$type %in% c("conventional Dendritic cells 2", "conventional Dendritic cells 1"),]$frac)
  colnames(maynard_census_fractions) <- c("type","mRNA")
  maynard_reads_fractions <- final[,c("type","mRNA")]
  maynard_census_fractions[maynard_census_fractions$type=="conventional Dendritic cells",]$mRNA<-mean(maynard_census_fractions[maynard_census_fractions$type %in% c("conventional Dendritic cells 2", "conventional Dendritic cells 1"),]$mRNA)
  
  
  rm(mat)
  
  return(list(maynard_census_fractions=maynard_census_fractions,
              maynard_reads_fractions=maynard_reads_fractions))
  
}


load_quantiseq <- function(){
  quantiseq <- fread("/home/Data/quanitseq.tsv")
  colnames(quantiseq) <- c("type","mRNA")
  quantiseq$type <- c("B cells","Macrophages","MacrophagesM2","Monocytes","Neutrophils","NK cells","T cells CD4","T cells CD8","T regulatory cells","Dendritic cells")
  
  return(quantiseq)
}

load_epic <- function(){
  load("/home/Data/epic.rda")
  epic <- data.frame(mRNA_cell_default)
  epic$type <- c("B cells","Macrophages","Monocytes","Neutrophils","NK cells" ,"T cells","T cells CD4","T cells CD8","T helper cells","T regulatory cells","otherCells","default")
  colnames(epic)<-c("mRNA","type")
  
  return(epic)
}

load_miltenyi <- function(){
  miltenyi <- data.frame(type=c("B cells", "Hematopoietic progenitor cells", "Monocytes", "conventional Dendritic cells", "plasmacytoid Dendritic cells","NK cells", "NKT cells", "T cells", "T cells CD4", "T cells CD8", "T regulatory cells"), mRNA=c(0.4, 0.5, 1, 2.5, 0.7, 0.3, 0.6, 1, 0.6, 0.4, 0.3))
  
  return(miltenyi)
}

# https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-6678/downloads?colourBy=metadata&metadata=inferred_cell_type_-_authors_labels
load_vento_tormo <- function(){
  meta_runs <- fread("/home/Data/Vento-Tormo/filereport_read_run_PRJEB25794_tsv.txt")[,c("library_name","read_count","run_accession")]
  meta_cells <- fread("/home/Data/Vento-Tormo/meta_ss2.txt")[,c("V1","annotation","location")]
  vento_census <- fread("/home/Data/Vento-Tormo/census_vento_tormo.tsv")
  meta_cells$ID <- gsub("#", "_", meta_cells$V1)
  meta_runs$library_name <- gsub("_p", "", meta_runs$library_name)
  meta <- merge(meta_cells,meta_runs, by.x="ID", by.y = "library_name")
  meta <- merge(meta, vento_census, by.x="run_accession", by.y="run")
  colnames(meta) <- c("id","id2","id3","cell_type","location","n_reads","run_accession","census")
  meta[,cell_count:=.N, by="cell_type"]
  ordered_types_blood <- c("conventional Dendritic cells 1", "conventional Dendritic cells 2", "Granulocytes", "Monocytes" ,"NK cells CD16 -","NK cells CD56bright CD16+/-","Plasma cells" ,"T cells","villous cytotrophoblast","decidual Macrophages 1"        
                     ,"NK cells proliferating","decidual NK cells 2","decidual NK cells 3")
  ordered_types_decidua <- c("conventional Dendritic cells 1","extravillous trophoblast", "Endothelial cells (f)" ,"Endothelial cells (m)", "Endothelial cells L","epithelial glandular cells 1","epithelial glandular cells 2",  "Granulocytes", "Hoffbauer cells","innate lymphocyte cells","NK cells CD16 -","NK cells CD56bright CD16+/-","Plasma cells" ,"syncytiotrophoblast","T cells","villous cytotrophoblast","decidual Macrophages 1"        
                           ,"decidual Macrophages 2","decidual Macrophages 3","NK cells proliferating","decidual NK cells 1","decidual NK cells 2","decidual NK cells 3","dP1" ,"dP2","stromal cells 1","stromal cells 2","stromal cells 3","Fibroblasts 1","Fibroblasts 2")
  ordered_types <- c("conventional Dendritic cells 1","conventional Dendritic cells 2","extravillous trophoblast", "Endothelial cells (f)" ,"Endothelial cells (m)", "Endothelial cells L","epithelial glandular cells 1","epithelial glandular cells 2",  "Granulocytes","Monocytes", "Hoffbauer cells","innate lymphocyte cells","NK cells CD16 -","NK cells CD56bright CD16+/-","Plasma cells" ,"syncytiotrophoblast","T cells","villous cytotrophoblast","decidual Macrophages 1"        
                     ,"decidual Macrophages 2","decidual Macrophages 3","NK cells proliferating","decidual NK cells 1","decidual NK cells 2","decidual NK cells 3","dP1" ,"dP2","stromal cells 1","stromal cells 2","stromal cells 3","Fibroblasts 1","Fibroblasts 2")
  
  
  vento_tormo_reads <- meta[, mean(n_reads), by=cell_type]
  colnames(vento_tormo_reads) <- c("type","mRNA")
  vento_tormo_reads <- vento_tormo_reads[order(type)]
  vento_tormo_reads$type <- ordered_types

  vento_tormo_census <- meta[, mean(census), by=cell_type]
  colnames(vento_tormo_census) <- c("type","mRNA")
  vento_tormo_census <- vento_tormo_census[order(type)]
  vento_tormo_census$type <- ordered_types
  
  
  # meta_blood <- meta[meta$location=="Blood"]
  # meta_decidua <- meta[meta$location=="Decidua"]
  
  #separate calculations of mean for both tissues
  # vento_tormo_blood_reads <- meta_blood[, mean(n_reads), by=cell_type]
  # colnames(vento_tormo_blood_reads) <- c("type","mRNA")
  # vento_tormo_blood_reads <- vento_tormo_blood_reads[order(type)]
  # vento_tormo_blood_reads$type <- ordered_types_blood
  # 
  # vento_tormo_blood_census <- meta_blood[, mean(census), by=cell_type]
  # colnames(vento_tormo_blood_census) <- c("type","mRNA")
  # vento_tormo_blood_census <- vento_tormo_blood_census[order(type)]
  # vento_tormo_blood_census$type <- ordered_types_blood
  # 
  # vento_tormo_decidua_reads <- meta_decidua[, mean(n_reads), by=cell_type]
  # colnames(vento_tormo_decidua_reads) <- c("type","mRNA")
  # vento_tormo_decidua_reads <- vento_tormo_decidua_reads[order(type)]
  # vento_tormo_decidua_reads$type <- ordered_types_decidua
  # 
  # vento_tormo_decidua_census <- meta_decidua[, mean(census), by=cell_type]
  # colnames(vento_tormo_decidua_census) <- c("type","mRNA")
  # vento_tormo_decidua_census <- vento_tormo_decidua_census[order(type)]
  # vento_tormo_decidua_census$type <- ordered_types_decidua
  # 
  # return(list(vento_tormo_blood_census=vento_tormo_blood_census,
  #             vento_tormo_blood_reads=vento_tormo_blood_reads,
  #             vento_tormo_decidua_census=vento_tormo_decidua_census,
  #             vento_tormo_decidua_reads=vento_tormo_decidua_reads))
  
  return(list(vento_tormo_census=vento_tormo_census,
              vento_tormo_reads=vento_tormo_reads))
}

load_vento_tormo_gene <- function(query_gene="ENSG00000126067"){
  genes <- fread("/home/Data/Vento-Tormo/E-MTAB-6678.expression_tpm.mtx_rows", header=F)[[1]]
  ordered_types <- c("conventional Dendritic cells 1","conventional Dendritic cells 2","extravillous trophoblast", "Endothelial cells (f)" ,"Endothelial cells (m)", "Endothelial cells L","epithelial glandular cells 1","epithelial glandular cells 2",  "Granulocytes","Monocytes", "Hoffbauer cells","innate lymphocyte cells","NK cells CD16 -","NK cells CD56bright CD16+/-","Plasma cells" ,"syncytiotrophoblast","T cells","villous cytotrophoblast","decidual Macrophages 1"        
                     ,"decidual Macrophages 2","decidual Macrophages 3","NK cells proliferating","decidual NK cells 1","decidual NK cells 2","decidual NK cells 3","dP1" ,"dP2","stromal cells 1","stromal cells 2","stromal cells 3","Fibroblasts 1","Fibroblasts 2")
  if(query_gene %in% genes){
    mat <- readMM("/home/Data/Vento-Tormo/E-MTAB-6678.expression_tpm.mtx")
    cells <- fread("/home/Data/Vento-Tormo/E-MTAB-6678.expression_tpm.mtx_cols", header=F)[[1]]
    dimnames(mat)<-list(genes, cells)
    
    query_expr <- data.frame(mat[which(query_gene == genes),])
    query_expr$run_accession <- rownames(query_expr)
    rm(mat)
    
    meta_runs <- fread("/home/Data/Vento-Tormo/filereport_read_run_PRJEB25794_tsv.txt")[,c("library_name","read_count","run_accession")]
    meta_cells <- fread("/home/Data/Vento-Tormo/meta_ss2.txt")[,c("V1","annotation","location")]
    meta_cells$ID <- gsub("#", "_", meta_cells$V1)
    meta_runs$library_name <- gsub("_p", "", meta_runs$library_name)
    meta <- merge(meta_cells,meta_runs, by.x="ID", by.y = "library_name")
    
    query_annotation <- merge(meta, query_expr, by="run_accession")
    colnames(query_annotation) <- c("run_accession","id1","id2","cell_type","location","n_reads","expression")
    
    expression_by_type <- query_annotation[,mean(expression), by=cell_type]
    colnames(expression_by_type) <- c("type","mRNA")
    expression_by_type <- expression_by_type[order(type)]
    expression_by_type$type <- ordered_types
    
    return(vento_tormo_query_gene=expression_by_type)
  }
}

load_monaco <- function(){
  cell_types <- c("Monocytes", "NK cells", "T cells CD8 central memory", "T cells CD4","T cells CD8","B Naive" ,"T cells CD4 central memory","MAIT" ,"T gd Vd2" ,"Neutrophils" ,"T gd non-Vd2" ,"Basophils LD", "Monocytes NC+I" ,"B cells"   ,"conventional Dendritic cells"   ,"plasmacytoid Dendritic cells" ,"Plasma cells") 
  fact <- c(22824.322779,  21456.905980,  17827.199342,  14262.064545,  10660.952413,  24162.648934,  15682.087676,  13365.625992,  28128.898864,  9546.742675,  25008.840273,  8610.360645,  22344.663722,  20837.569855,  57322.183941,  28488.194770,  325800.989149)
  monaco <- data.frame(type=cell_types, mRNA=fact)
  
  return(monaco)
}


#The counts were first re-scaled by library-size per cell. As the expression-cutoff i chose 0, since we are not yet sure how this factor works exactly and also I believe it is chosen on a TPM distribution and not relative expression counts (there are almost never counts with log10(count) > 0.1; 0.1=default in monocle for TPM/FPKM).
load_hao <- function(){
  df <- fread("/home/Data/Hao (CITEseq-PBMC)/census_meta.tsv")
  hao_census <-df[,mean(census_counts),by=celltype.l2]
  colnames(hao_census) <- c("type","mRNA")
  hao_census <- hao_census[order(type)]
  ordered_types <-
    c(
      "ASDC",
      "B cells intermediate",
      "B cells memory",
      "B cells",
      "Monocytes CD14",
      "Monocytes CD16",
      "T cells CD4 cytotoxic activity",
      "T cells CD4",
      "T cells CD4 proliferating",
      "T cells CD4 central memory",
      "T cells CD4 effector memory",
      "T cells CD8",
      "T cells CD8 proliferating",
      "T cells CD8 central memory",
      "T cells CD8 effector memory",
      "Doublet",
      "Eryth",
      "HPSC",
      "innate lymphocyte cells",
      "MAI T cells",
      "NK cells",
      "NK cells proliferating",
      "NK cells CD56bright CD16+/-",
      "Plasma cells",
      "Platelet",
      "T regulatory cells",
      "conventional Dendritic cells 1",
      "conventional Dendritic cells 2",
      "T cells double negative",
      "T cells gamma delta",
      "plasmacytoid Dendritic cells"
    )
  
  hao_reads <- df[,mean(nCount_RNA), by=celltype.l2]
  colnames(hao_reads) <- c("type","mRNA")
  hao_reads <- hao_reads[order(type)]
  
  hao_reads$type <- ordered_types
  hao_census$type <- ordered_types
  
  return(list(hao_reads=hao_reads,
              hao_census=hao_census))
}

load_hao_psmb2 <- function(){
  mean_expr <- fread("/home/Data/Hao (CITEseq-PBMC)/mean_expr_psmb2.tsv")
  mean_expr$V1 <- NULL
  colnames(mean_expr) <- c("type","mRNA")
  mean_expr <- mean_expr[order(type)]
  ordered_types <-
    c(
      "ASDC",
      "B cells intermediate",
      "B cells memory",
      "B cells",
      "Monocytes CD14",
      "Monocytes CD16",
      "T cells CD4 cytotoxic activity",
      "T cells CD4",
      "T cells CD4 proliferating",
      "T cells CD4 central memory",
      "T cells CD4 effector memory",
      "T cells CD8",
      "T cells CD8 proliferating",
      "T cells CD8 central memory",
      "T cells CD8 effector memory",
      "Doublet",
      "Eryth",
      "HPSC",
      "innate lymphocyte cells",
      "MAI T cells",
      "NK cells",
      "NK cells proliferating",
      "NK cells CD56bright CD16+/-",
      "Plasma cells",
      "Platelet",
      "T regulatory cells",
      "conventional Dendritic cells 1",
      "conventional Dendritic cells 2",
      "T cells double negative",
      "T cells gamma delta",
      "plasmacytoid Dendritic cells"
    )
  mean_expr$type <- ordered_types
  return(hao_gene=mean_expr)
}


maynard <- load_maynard()
maynard_reads <- maynard$maynard_reads_fractions
maynard_census <- maynard$maynard_census_fractions

vento <- load_vento_tormo()
# vento_blood_reads <- vento$vento_tormo_blood_reads
# vento_blood_census <- vento$vento_tormo_blood_census
# vento_decidua_reads <- vento$vento_tormo_decidua_reads
# vento_decidua_census <- vento$vento_tormo_decidua_census

vento_reads <- vento$vento_tormo_reads
vento_census <- vento$vento_tormo_census

vento_psmb2 <- load_vento_tormo_gene()

quantiseq <- load_quantiseq()
epic <- load_epic()
monaco <- load_monaco()
miltenyi <- load_miltenyi()

hao <- load_hao()
hao_reads <- hao$hao_reads
hao_census <- hao$hao_census

hao_psmb2 <- load_hao_psmb2()



