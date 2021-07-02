library("quantiseqr")
library("dplyr")
library("ggplot2")
library("data.table")
library("tidyr")
library("tibble")
library("GEOquery")
library("reshape2")
#library("SummarizedExperiment")
library("Matrix")
library("Seurat")
library("ggpubr")


##### quantiseq vignette #####

tpminfo_GSE107572 <- getGEOSuppFiles("GSE107572",
                                     baseDir = tempdir(),
                                     filter_regex = "GSE107572_tpm_PBMC_RNAseq"
)
tpm_location <- rownames(tpminfo_GSE107572)[1]
tpmdata <- read.table(tpm_location, header = TRUE)

tpm_genesymbols <- tpmdata$GENE
tpmdata <- as.matrix(tpmdata[, -1])
rownames(tpmdata) <- tpm_genesymbols

ti_PBMCs <- quantiseqr::run_quantiseq(
  expression_data = tpmdata,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = FALSE,
  scale_mRNA = TRUE
)


GEOid <- "GSE107572"
gds <- getGEO(GEOid)

GEOinfo <- pData(gds[[1]])
FACSdata <- data.frame(
  B.cells = GEOinfo$`b cells:ch1`,
  T.cells.CD4 = GEOinfo$`cd4+ t cells:ch1`,
  T.cells.CD8 = GEOinfo$`cd8+ t cells:ch1`,
  Monocytes = GEOinfo$`monocytes:ch1`,
  Dendritic.cells = GEOinfo$`myeloid dendritic cells:ch1`,
  NK.cells = GEOinfo$`natural killer cells:ch1`,
  Neutrophils = GEOinfo$`neutrophils:ch1`,
  Tregs = GEOinfo$`tregs:ch1`
)
FACSdata$Sample <- gsub(
  "Blood-derived immune-cell mixture from donor ", "pbmc", GEOinfo$title
)

ti_PBMCs$Sample <- gsub("_.*$", "", sub("_", "", ti_PBMCs$Sample))


ti <- gather(ti_PBMCs, celltype, fraction, B.cells:Other)
#ti$classification <- "quantiseq"

facs <- gather(FACSdata, celltype, fraction, B.cells:Tregs)
#facs$classification <- "FACS"
facs$fraction <- as.numeric(facs$fraction)

combined_data_fino <- merge(ti, facs, by = c("Sample","celltype"))
colnames(combined_data_fino) <- c("Sample","celltype","fraction_quantiseq","fraction_facs")
#combined_data_fino$dataset <- "Finotello"


##### Hoek dataset #####

hoek_tpm <- read.table("/home/Data/Hoek/Hoek_data/Hoek_data/Hoek_PBMC_TPM.txt", header=T, row.names = 1)

hoek_quanti_orig <- quantiseqr::run_quantiseq(
  expression_data = hoek_tpm,
  signature_matrix = "TIL10",
  is_arraydata = FALSE,
  is_tumordata = FALSE,
  scale_mRNA = TRUE
)
hoek_quanti <- gather(hoek_quanti_orig, celltype, fraction, B.cells:Other)


hoek_facs_orig <- read.table("/home/Data/Hoek/Hoek_data/Hoek_data/Hoek_PBMC_FACS.txt", header=T)
colnames(hoek_facs_orig) <- c("Sample","T.cells.CD4","Monocytes","B.cells","Dendritic.cells","NK.cells")
hoek_facs <- gather(hoek_facs_orig, celltype, fraction, T.cells.CD4:NK.cells)


combined_data_hoek <- merge(hoek_quanti, hoek_facs, by = c("Sample","celltype"))
colnames(combined_data_hoek) <- c("Sample","celltype","fraction_quantiseq","fraction_facs")
#combined_data_hoek$dataset <- "Hoek"


##### plot celltype fractions of two datasets #####

p1<-ggscatter(combined_data_hoek, x="fraction_quantiseq",y="fraction_facs", 
             cor.coef = T, facet.by = "celltype", color = "celltype",
             scales="free", palette = "jco", conf.int = T, size=4.5, 
             title="Celltype fractions for Hoek PBMC dataset (quantiseq vs FACS)")
p1<-p1+geom_smooth(method = "lm", se=F, color="black")+geom_abline(color="lightgrey",linetype="dashed")
ggsave("plots/Hoek_fractions.png",plot = p1, scale = 2)

p2<-ggscatter(combined_data_fino, x="fraction_quantiseq",y="fraction_facs", 
             cor.coef = T, facet.by = "celltype", color = "celltype",
             scales="free", palette = "jco", conf.int = T, size=4.5, 
             title="Celltype fractions for Finotello PBMC dataset (quantiseq vs FACS)")
p2<-p2+geom_smooth(method = "lm", se=F, color="black")+geom_abline(color="lightgrey",linetype="dashed")
ggsave("plots/Finotello_fractions.png",plot = p2, scale = 2)

combined_data_hoek$dataset <- "Hoek"
combined_data_fino$dataset <- "Finotello"
data <- rbind(combined_data_hoek, combined_data_fino)

p3<-ggscatter(data, x="fraction_quantiseq",y="fraction_facs", shape="dataset",
            cor.coef = T, facet.by = "celltype", color = "celltype",
            scales="free", palette = "jco", conf.int = T, size=4.5, 
            title="Celltype fractions for two PBMC datasets (quantiseq vs FACS)")
p3<-p3+geom_smooth(method = "lm", se=F, color="black")+geom_abline(color="lightgrey",linetype="dashed")
ggsave("plots/Finotello_Hoek_fractions.png",plot = p3, scale = 2)


##### work with immunedeconv results #####

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# gaussian kernel densitiy estimation (https://stackoverflow.com/questions/64235786/gaussian-kernel-density-estimation-in-r)
GK <- function(u) {(1/sqrt(2*pi))*exp(-(u^2)/2)}

deconv <- fread("../Data/Finotello/results_deconvolution.tsv")

ggplot(deconv, aes(x=value, y=sample,fill=cell_type))+
  geom_bar(stat="identity")+
  facet_grid(~method, scale="free")

##### Hao dataset #####

#HTO5p.mtx<-read.csv(gzfile("../Data/Hao (CITEseq-PBMC)/GSM5008742_HTO_5P-matrix.mtx.gz"), sep=" ")
#HTO5p.feat<-read.csv(gzfile("../Data/Hao (CITEseq-PBMC)/GSM5008742_HTO_5P-features.tsv.gz"), sep=" ")
#HTO5P.barc<-read.csv(gzfile("../Data/Hao (CITEseq-PBMC)/GSM5008742_HTO_5P-barcodes.tsv.gz"), sep=" ")

source("scripts/census.R")
matrix_dir = "/home/Data/Hao (CITEseq-PBMC)/hto_5p/"
mat<-Read10X(matrix_dir)
out_paper <- census(matrix = mat, ncores = 1, method = "paper")



