source("/home/simulator/scripts/census.R")

load_trivaglini <- function(){
  meta_file <- "/home/Data/Travaglini/obs_extended.csv"
  
  meta <- fread(meta_file)
  meta$free_annotation[meta$free_annotation=="Adventitial Fibroblast_P1"] <- "Fibroblast"
  meta$free_annotation[meta$free_annotation=="Adventitial Fibroblast_P2"] <- "Fibroblast"
  meta$free_annotation[meta$free_annotation=="Adventitial Fibroblast_P3"] <- "Fibroblast"
  meta$free_annotation[meta$free_annotation=="Airway Smooth Muscle_P1"] <- "Smooth muscle cell Airway"
  meta$free_annotation[meta$free_annotation=="Airway Smooth Muscle_P2"] <- "Smooth muscle cell Airway"
  meta$free_annotation[meta$free_annotation=="Airway Smooth Muscle_P3"] <- "Smooth muscle cell Airway"
  meta$free_annotation[meta$free_annotation=="Alveolar Epithelial Type 1_P1"] <- "Epithelial cell Alveolar 1"
  meta$free_annotation[meta$free_annotation=="Alveolar Epithelial Type 1_P2"] <- "Epithelial cell Alveolar 1"
  meta$free_annotation[meta$free_annotation=="Alveolar Epithelial Type 1_P3"] <- "Epithelial cell Alveolar 1"
  meta$free_annotation[meta$free_annotation=="Alveolar Epithelial Type 2_P1"] <- "Epithelial cell Alveolar 2"
  meta$free_annotation[meta$free_annotation=="Alveolar Epithelial Type 2_P2"] <- "Epithelial cell Alveolar 2"
  meta$free_annotation[meta$free_annotation=="Alveolar Epithelial Type 2_P3"] <- "Epithelial cell Alveolar 2"
  meta$free_annotation[meta$free_annotation=="Alveolar Fibroblast_P1"] <- "Fibroblast"
  meta$free_annotation[meta$free_annotation=="Alveolar Fibroblast_P2"] <- "Fibroblast"
  meta$free_annotation[meta$free_annotation=="Alveolar Fibroblast_P3"] <- "Fibroblast"
  meta$free_annotation[meta$free_annotation=="Artery_P1"] <- "Artery 1"
  meta$free_annotation[meta$free_annotation=="Artery_P2"] <- "Artery 2"
  meta$free_annotation[meta$free_annotation=="Artery_P3"] <- "Artery 3"
  meta$free_annotation[meta$free_annotation=="B_P1"] <- "B cells"
  meta$free_annotation[meta$free_annotation=="B_P2"] <- "B cells"
  meta$free_annotation[meta$free_annotation=="B_P3"] <- "B cells"
  meta$free_annotation[meta$free_annotation=="Basal_P1"] <- "Basal"
  meta$free_annotation[meta$free_annotation=="Basal_P2"] <- "Basal"
  meta$free_annotation[meta$free_annotation=="Basal_P3"] <- "Basal"
  meta$free_annotation[meta$free_annotation=="Basophil/Mast 1_P1"] <- "Mast cells"
  meta$free_annotation[meta$free_annotation=="Basophil/Mast 1_P2"] <- "Mast cells"
  meta$free_annotation[meta$free_annotation=="Basophil/Mast 1_P3"] <- "Mast cells"
  meta$free_annotation[meta$free_annotation=="Bronchial Vessel 1_P1"] <- "Bronchial Vessel"
  meta$free_annotation[meta$free_annotation=="Capillary Aerocyte_P1"] <- "Capillary Aerocyte 1"
  meta$free_annotation[meta$free_annotation=="Capillary Aerocyte_P2"] <- "Capillary Aerocyte 2"
  meta$free_annotation[meta$free_annotation=="Capillary Aerocyte_P3"] <- "Capillary Aerocyte 3"
  meta$free_annotation[meta$free_annotation=="Capillary Intermediate 1_P2"] <- "Capillary Intermediate"
  meta$free_annotation[meta$free_annotation=="Capillary_P1"] <- "Capillary 1"
  meta$free_annotation[meta$free_annotation=="Capillary_P2"] <- "Capillary 2"
  meta$free_annotation[meta$free_annotation=="Capillary_P3"] <- "Capillary 3"
  meta$free_annotation[meta$free_annotation=="CD4+ Memory/Effector T_P1"] <- "T cells CD4 central memory"
  meta$free_annotation[meta$free_annotation=="CD4+ Naive T_P1"] <- "T cells CD4"
  meta$free_annotation[meta$free_annotation=="CD4+ Naive T_P2"] <- "T cells CD4"
  meta$free_annotation[meta$free_annotation=="CD8+ Memory/Effector T_P1"] <- "T cells CD8 central memory"
  meta$free_annotation[meta$free_annotation=="CD8+ Naive T_P1"] <- "T cells CD8"
  meta$free_annotation[meta$free_annotation=="CD8+ Naive T_P2"] <- "T cells CD8"
  meta$free_annotation[meta$free_annotation=="Ciliated_P1"] <- "Ciliated"
  meta$free_annotation[meta$free_annotation=="Ciliated_P2"] <- "Ciliated"
  meta$free_annotation[meta$free_annotation=="Ciliated_P3"] <- "Ciliated"
  meta$free_annotation[meta$free_annotation=="Classical Monocyte_P1"] <- "Monocytes"
  meta$free_annotation[meta$free_annotation=="Club_P1"] <- "Club cells"
  meta$free_annotation[meta$free_annotation=="Club_P2"] <- "Club cells"
  meta$free_annotation[meta$free_annotation=="Club_P3"] <- "Club cells"
  meta$free_annotation[meta$free_annotation=="Dendritic_P1"] <- "Dendritic cells"
  meta$free_annotation[meta$free_annotation=="Differentiating Basal_P3"] <- "Differentiating Basal"
  meta$free_annotation[meta$free_annotation=="Fibromyocyte_P3"] <- "Fibromyocyte"
  meta$free_annotation[meta$free_annotation=="Goblet_P1"] <- "Goblet"
  meta$free_annotation[meta$free_annotation=="Goblet_P2"] <- "Goblet"
  meta$free_annotation[meta$free_annotation=="Goblet_P3"] <- "Goblet"
  meta$free_annotation[meta$free_annotation=="IGSF21+ Dendritic_P2"] <- "Dendritic cells IGSF21+ 1"
  meta$free_annotation[meta$free_annotation=="IGSF21+ Dendritic_P3"] <- "Dendritic cells IGSF21+ 2"
  meta$free_annotation[meta$free_annotation=="Intermediate Monocyte_P2"] <- "Monocyte intermediate 1"
  meta$free_annotation[meta$free_annotation=="Intermediate Monocyte_P3"] <- "Monocyte intermediate 2"
  meta$free_annotation[meta$free_annotation=="Ionocyte_P3"] <- "Ionocyte"
  meta$free_annotation[meta$free_annotation=="Lipofibroblast_P1"] <- "Lipofibroblast"
  meta$free_annotation[meta$free_annotation=="Lymphatic_P1"] <- "Lymphatic 1"
  meta$free_annotation[meta$free_annotation=="Lymphatic_P2"] <- "Lymphatic 2"
  meta$free_annotation[meta$free_annotation=="Lymphatic_P3"] <- "Lymphatic 3"
  meta$free_annotation[meta$free_annotation=="Macrophage_P2"] <- "Macrophages"
  meta$free_annotation[meta$free_annotation=="Macrophage_P3"] <- "Macrophages"
  meta$free_annotation[meta$free_annotation=="Myeloid Dendritic Type 2_P3"] <- "Dendritic cells conventional"
  meta$free_annotation[meta$free_annotation=="Myofibroblast_P2"] <- "Myofibroblast 1"
  meta$free_annotation[meta$free_annotation=="Myofibroblast_P3"] <- "Myofibroblast 2"
  meta$free_annotation[meta$free_annotation=="Natural Killer T_P2"] <- "NK cells T"
  meta$free_annotation[meta$free_annotation=="Natural Killer T_P3"] <- "NK cells T"
  meta$free_annotation[meta$free_annotation=="Natural Killer_P1"] <- "NK cells"
  meta$free_annotation[meta$free_annotation=="Natural Killer_P2"] <- "NK cells"
  meta$free_annotation[meta$free_annotation=="Natural Killer_P3"] <- "NK cells"
  meta$free_annotation[meta$free_annotation=="Neuroendocrine_P1"] <- "Neuroendocrine 1"
  meta$free_annotation[meta$free_annotation=="Neuroendocrine_P3"] <- "Neuroendocrine 2"
  meta$free_annotation[meta$free_annotation=="Neutrophil_P1"] <- "Neutrophils"
  meta$free_annotation[meta$free_annotation=="Neutrophil_P2"] <- "Neutrophils"
  meta$free_annotation[meta$free_annotation=="Neutrophil_P3"] <- "Neutrophils"
  meta$free_annotation[meta$free_annotation=="Nonclassical Monocyte_P1"] <- "Monocyte non-conventional"
  meta$free_annotation[meta$free_annotation=="Nonclassical Monocyte_P2"] <- "Monocyte non-conventional"
  meta$free_annotation[meta$free_annotation=="Pericyte_P1"] <- "Pericyte"
  meta$free_annotation[meta$free_annotation=="Pericyte_P2"] <- "Pericyte"
  meta$free_annotation[meta$free_annotation=="Pericyte_P3"] <- "Pericyte"
  meta$free_annotation[meta$free_annotation=="Plasma_P3"] <- "Plasma cells"
  meta$free_annotation[meta$free_annotation=="Plasmacytoid Dendritic_P1"] <- "Dendritic cells plasmacytoid"
  meta$free_annotation[meta$free_annotation=="Plasmacytoid Dendritic_P2"] <- "Dendritic cells plasmacytoid"
  meta$free_annotation[meta$free_annotation=="Plasmacytoid Dendritic_P3"] <- "Dendritic cells plasmacytoid"
  meta$free_annotation[meta$free_annotation=="Proliferating NK/T_P2"] <- "NK cells proliferating"
  meta$free_annotation[meta$free_annotation=="Proliferating NK/T_P3"] <- "NK cells proliferating"
  meta$free_annotation[meta$free_annotation=="Signaling Alveolar Epithelial Type 2_P1"] <- "Epithelial cell singaling Alveolar"
  meta$free_annotation[meta$free_annotation=="Signaling Alveolar Epithelial Type 2_P3"] <- "Epithelial cell singaling Alveolar"
  meta$free_annotation[meta$free_annotation=="Vascular Smooth Muscle_P1"] <- "Smooth muscle cell Muscle"
  meta$free_annotation[meta$free_annotation=="Vascular Smooth Muscle_P2"] <- "Smooth muscle cell Muscle"
  meta$free_annotation[meta$free_annotation=="Vascular Smooth Muscle_P3"] <- "Smooth muscle cell Muscle"
  meta$free_annotation[meta$free_annotation=="Vein_P2"] <- "Vein"
  
  trivaglini_reads <- meta[,mean(nReads),by="free_annotation"]
  colnames(trivaglini_reads) <-c("type","mRNA")
  trivaglini_reads <- rbind(trivaglini_reads, list("T cells", mean(trivaglini_reads[trivaglini_reads$type %in% c("T cells CD4 central memory","T cells CD8 central memory","T cells CD4","T cells CD8"),]$mRNA)))
  
  trivaglini_spike <- meta[,mean(percent.ercc),by="free_annotation"]
  colnames(trivaglini_spike) <-c("type","mRNA")
  trivaglini_spike <- rbind(trivaglini_spike, list("T cells", mean(trivaglini_spike[trivaglini_spike$type %in% c("T cells CD4 central memory","T cells CD8 central memory","T cells CD4","T cells CD8"),]$mRNA)))
  
  ercc_cols <- grep("ERCC-",colnames(meta))
  meta_red <- data.frame(meta)[ercc_cols]
  meta_red$ercc_sum <- apply(meta_red,1,sum)
  meta_red$type <-meta$free_annotation
  meta_red$total_reads <- meta$nReads
  meta_red$genes <- meta$nGene
  meta_red <- as.data.table(meta_red)
  # i am adding a pseudocount of +1 here, to not have 0 in the Zähler for most of the cell-types
  meta_scaling <- meta_red[,lapply(.SD, function(x){mean(genes/na.omit(x+1))}),.SDcols=grep("ERCC",colnames(meta_red)), by="type"]
  
  # remove spike-ins with 0 counts
  meta_scaling <- as.data.frame(meta_scaling)[,colSums(meta_scaling!=0)>0]
  
  trivaglini_ercc_total <- meta_red[,mean(ercc_sum),by="type"]
  colnames(trivaglini_ercc_total) <-c("type","mRNA")
  trivaglini_ercc_total <- rbind(trivaglini_ercc_total, list("T cells", mean(trivaglini_ercc_total[trivaglini_ercc_total$type %in% c("T cells CD4 central memory","T cells CD8 central memory","T cells CD4","T cells CD8"),]$mRNA)))
  
  
  trivaglini_genes_per_ercc <- meta_red[,mean(genes/ercc_sum), by="type"]
  colnames(trivaglini_genes_per_ercc) <-c("type","mRNA")
  trivaglini_genes_per_ercc <- rbind(trivaglini_genes_per_ercc, list("T cells", mean(trivaglini_genes_per_ercc[trivaglini_genes_per_ercc$type %in% c("T cells CD4 central memory","T cells CD8 central memory","T cells CD4","T cells CD8"),]$mRNA)))
  
  
  return(list(trivaglini_reads=trivaglini_reads,
              trivaglini_spike=trivaglini_spike,
              trivaglini_ercc=meta_scaling,
              trivaglini_ercc_total=trivaglini_ercc_total,
              trivaglini_genes_per_ercc=trivaglini_genes_per_ercc))
}

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
  colnames(maynard_census_fractions) <- c("type","mRNA")
  maynard_census_fractions <- rbind(maynard_census_fractions, list("convetional Dendritic cells", mean(maynard_census_fractions[maynard_census_fractions$type %in% c("conventional Dendritic cells 1", "conventional Dendritic cells 2"),]$mRNA)))
  maynard_census_fractions <- rbind(maynard_census_fractions, list("Fibroblasts", mean(maynard_census_fractions[maynard_census_fractions$type %in% c("Fibroblast adventitial", "Fibroblast alevolar"),]$mRNA)))
  
  
  maynard_reads_fractions <- final[,c("type","mRNA")]
  maynard_reads_fractions <- rbind(maynard_reads_fractions, list("convetional Dendritic cells", mean(maynard_reads_fractions[maynard_reads_fractions$type %in% c("conventional Dendritic cells 2", "conventional Dendritic cells 1"),]$mRNA)))
  maynard_reads_fractions <- rbind(maynard_reads_fractions, list("Fibroblasts", mean(maynard_reads_fractions[maynard_reads_fractions$type %in% c("Fibroblast adventitial", "Fibroblast alevolar"),]$mRNA)))
  
  
  rm(mat)
  
  return(list(maynard_census_fractions=maynard_census_fractions,
              maynard_reads_fractions=maynard_reads_fractions))
  
}


load_quantiseq <- function(){
  quantiseq <- fread("/home/Data/quanitseq.tsv")
  colnames(quantiseq) <- c("type","mRNA")
  quantiseq$type <- c("B cells","Macrophages","MacrophagesM2","Monocytes","Neutrophils","NK cells","T cells CD4","T cells CD8","T regulatory cells","Dendritic cells")
  
  quantiseq <- rbind(quantiseq, list("T cells", mean(quantiseq[quantiseq$type %in% c("T cells CD4", "T cells CD8", "T regulatory cells"),]$mRNA)))

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
  vento_tormo_reads<-rbind(vento_tormo_reads, list("NK cells", mean(vento_tormo_reads[vento_tormo_reads$type %in% c("NK cells CD16 -", "NK cells CD56bright CD16+/-","NK cells proliferating"),]$mRNA)))
  vento_tormo_reads<-rbind(vento_tormo_reads, list("Fibroblasts", mean(vento_tormo_reads[vento_tormo_reads$type %in% c("Fibroblasts 1", "Fibroblasts 2"),]$mRNA)))
  vento_tormo_reads<-rbind(vento_tormo_reads, list("Dendritic cells", mean(vento_tormo_reads[vento_tormo_reads$type %in% c("conventional Dendritic cells 1", "conventional Dendritic cells 2"),]$mRNA)))
  
  vento_tormo_census <- meta[, mean(census), by=cell_type]
  colnames(vento_tormo_census) <- c("type","mRNA")
  vento_tormo_census <- vento_tormo_census[order(type)]
  vento_tormo_census$type <- ordered_types
  vento_tormo_census<-rbind(vento_tormo_census, list("NK cells", mean(vento_tormo_census[vento_tormo_census$type %in% c("NK cells CD16 -", "NK cells CD56bright CD16+/-","NK cells proliferating"),]$mRNA)))
  vento_tormo_census<-rbind(vento_tormo_census, list("Fibroblasts", mean(vento_tormo_census[vento_tormo_census$type %in% c("Fibroblasts 1", "Fibroblasts 2"),]$mRNA)))
  vento_tormo_census<-rbind(vento_tormo_census, list("Dendritic cells", mean(vento_tormo_census[vento_tormo_census$type %in% c("conventional Dendritic cells 1", "conventional Dendritic cells 2"),]$mRNA)))
  
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
    expression_by_type<-rbind(expression_by_type, list("NK cells", mean(expression_by_type[expression_by_type$type %in% c("NK cells CD16 -", "NK cells CD56bright CD16+/-","NK cells proliferating"),]$mRNA)))
    expression_by_type<-rbind(expression_by_type, list("Fibroblasts", mean(expression_by_type[expression_by_type$type %in% c("Fibroblasts 1", "Fibroblasts 2"),]$mRNA)))
    expression_by_type<-rbind(expression_by_type, list("Dendritic cells", mean(expression_by_type[expression_by_type$type %in% c("conventional Dendritic cells 1", "conventional Dendritic cells 2"),]$mRNA)))
    
    
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
  
  hao_reads <- rbind(hao_reads, list("Monocytes", mean(hao_reads[hao_reads$type %in% c("Monocytes CD14", "Monocytes CD16"),]$mRNA)))
  hao_census <- rbind(hao_census, list("Monocytes", mean(hao_census[hao_census$type %in% c("Monocytes CD14", "Monocytes CD16"),]$mRNA)))
  
    
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
  
  mean_expr <- rbind(mean_expr, list("Monocytes", mean(mean_expr[mean_expr$type %in% c("Monocytes CD14", "Monocytes CD16"),]$mRNA)))
  
  
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

trivaglini <- load_trivaglini()
trivaglini_reads <- trivaglini$trivaglini_reads
trivaglini_spike <- trivaglini$trivaglini_spike
trivaglini_ercc <- trivaglini$trivaglini_ercc
trivaglini_ercc_total <- trivaglini$trivaglini_ercc_total
trivaglini_genes_per_ercc <- trivaglini$trivaglini_genes_per_ercc

