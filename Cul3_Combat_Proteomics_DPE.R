
library(tidyverse)
library(qPLEXanalyzer)
library(UniProt.ws)
library(dplyr)
library(sva)
library(limma)



dir.create("./Proteomics_Output_Plots_Period_Region")

### DATA PROCESSING ###

# Reading input data sheet with TMT labelled peptide intensities
dat<-read.table("./Peptide_intensities.txt", header = T,  sep = "\t", fill = TRUE)

# Filter the proteins which have unique peptide sequence counts larger than 2
filter_protein <- dat %>%
  dplyr::select(PROTEIN, SEQUENCE) %>%
  distinct(PROTEIN, SEQUENCE) %>%
  count(PROTEIN) %>%
  filter(n >= 2) %>%
  pull(PROTEIN) %>%
  as.character()

dat_filt <- dat%>%
  filter(PROTEIN %in% filter_protein) %>%
  dplyr::select(-starts_with("m"))%>%
  dplyr::select(
    PROTEIN, SEQUENCE, UNIQUE,starts_with("norm")
  ) %>%
  mutate(PROTEIN = as.character(PROTEIN))

# Reading the meta data sheet where contains information for samples
# Mapping the spectral counts to TMT labels and samples

meta_P7<-readxl::read_xlsx("meta_data.xlsx", sheet = 1)
meta_P7<- meta_P7[order(meta_P7$SampleGroup),]
meta_P7<-meta_P7%>%
  mutate(SampleName = paste(SampleName, lapply(str_split(meta_P7$SampleID,"-"), function(x)x[1]),sep = "_"))

names(dat_filt)<-c("PROTEIN", "SEQUENCE", "UNIQUE", "126_C1", "127N_C2", "127C_C2", 
                   "128N_C2", "128C_C2", "129N_C4", "129C_C4", "130N_C1", "130C_pooled", "131N_pooled")

# Removing contaminants
dat_filt<-dat_filt%>%
  filter(! str_detect(PROTEIN, "contaminant|Reverse"))

# Reordering the data
dat_filt_match<-dat_filt[,4:13]
dat_filt_match<-dat_filt_match[match(meta_P7$SampleName, colnames(dat_filt_match))]
dat_sort<-cbind(dat_filt[,1:3], dat_filt_match)

### Proteomics Analysis ###

# Creating MS proteomics expression object
exp2=list(intensities=dat_sort, metadata=meta_P7)
Mnset<-convertToMSnset(exp2$intensities, metadata = exp2$metadata, indExpData = c(4:13), 
                       Sequences = 2, Accessions = 1)
exprs(Mnset)<-exprs(Mnset)+0.01
proteins <- unique(fData(Mnset)$Accessions)
columns <- c("ENTRY-NAME", "PROTEIN-NAMES", "GENES")

# Annotation of the filtered proteins with uniprot
ms <- UniProt.ws::UniProt.ws(taxId = 10090)
mouse_anno <- UniProt.ws::select(ms, proteins, columns, "UNIPROTKB") %>%
  mutate(GeneSymbol = gsub(" .*", "", GENES)) %>%
  dplyr::select(
    Accessions = "UNIPROTKB", Gene = "ENTRY-NAME",
    Description = "PROTEIN-NAMES", GeneSymbol
  )
# Summarize peptides intensities to protein level
MSnset_Sum <- summarizeIntensities(Mnset, sum, mouse_anno)
unorm_sum_intensity<-exprs(MSnset_Sum)

data_start <- unorm_sum_intensity
data_start<-data_start[, -c(5:6)]

data_raw <- na.omit(data_start)
data_raw<- as.data.frame(data_raw)


# Sample Loading Normalization
exp1_raw <- data_raw[c(1:4)]
exp2_raw <- data_raw[c(5:8)]
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw)))
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")

data_sl <- cbind(exp1_sl, exp2_sl)

data_sl_pca<-data_sl


# Removing Batch Effect using combat
design <- meta_PCA
design$SampleGroup <- factor(design$SampleGroup, levels = c("Cul3 HET", "WT"))
mod <- model.matrix(~ SampleGroup + Gender, data = design)
batch <- design$Batch
data_sl=as.matrix(data_sl)
data_combat <- ComBat(dat = data_sl, batch = batch, mod = mod, par.prior = TRUE)
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x > 0)), ]
data_combat<-as.data.frame(data_combat)



# Differential protein expression analysis using limma 
mod_mat <- mod
colnames(mod_mat) <- make.names(colnames(mod_mat))
lm_fit <- lmFit(log2(data_combat), design = mod_mat)
contr_mat <- makeContrasts(
  HET = -SampleGroupWT,
  levels = colnames(mod_mat)
)
contr_fit <- contrasts.fit(lm_fit, contr_mat)
ebayes <- eBayes(contr_fit)
tt <- topTable(ebayes, coef = "HET", number = Inf)

mouse_anno_sort_limma=mouse_anno[match(rownames(tt), mouse_anno$Accessions),]
tt_limma<-cbind(tt, mouse_anno_sort_limma)
DPE_limma=tt_limma[tt_limma$adj.P.Val < 0.1,]
write.csv(tt_limma, "Limma_filter.csv")







