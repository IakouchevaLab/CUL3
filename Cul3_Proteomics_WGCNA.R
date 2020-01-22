library(tidyverse)
library(limma)
library(edgeR)
library(WGCNA)
library(plyr)
enableWGCNAThreads()
allowWGCNAThreads()


# Data Processing 
# Reading the proteomics input data after removing batch effect if it exists
# Reading the metadata sheet for corresponding proteomics input
Proteomics_input_data <- readxl::read_xlsx("Proteomics_Input_Data.xlsx", sheet = 1) %>%
as.data.frame(()
metadata_input <- readxl::read_xlsx("MetaData.xlsx", sheet = 1)
rownames(Proteomics_input_data) <- Proteomics_input_data$Protein_ID; Proteomics_input_data$Protein_ID=NULL


# calculate the soft threshold
powers = c(seq(1,9,by=1),seq(10,30,by=2))
suggested_sft = pickSoftThreshold(data= t(Proteomics_input_data), 
                            networkType = "signed", corFnc="bicor",verbose=2,
                            powerVector=powers,blockSize = 20000)
pdf("sft_cx_bicor_Prot_signed.pdf", height = 8, width = 8)
par(mfrow=c(1,2))
plot(suggested_sft$fitIndices[,1], 
     -sign(suggested_sft$fitIndices[,3])*suggested_sft$fitIndices[,2],
     xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(suggested_sft$fitIndices[,1], 
     -sign(suggested_sft$fitIndices[,3])*suggested_sft$fitIndices[,2], 
     labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", 
     ylab="Scale free R^2")
abline(h=0.8, col="black")
plot(suggested_sft$fitIndices[,1], suggested_sft$fitIndices[,5], 
     xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(suggested_sft$fitIndices[,1], suggested_sft$fitIndices[,5], 
     labels = powers, cex = 0.7, col="black")
dev.off()

# visualising the result by the dendrogram
pdf("Prot_sft_bicor_signed.pdf", height = 8, width = 8)
net = blockwiseModules(datExpr=t(Proteomics_input_data[, -9]), maxBlockSize=20000,networkType="signed",corType="bicor",  
                       power = 18, mergeCutHeight= 0.1, minModuleSize= 40, pamStage=TRUE, reassignThreshold=1e-6, 
                       saveTOMFileBase="Prot_WGCNA_signed_bicor_Prot", saveTOMs=TRUE, 
                       verbose = Inf, deepSplit=2)
dev.off()
save(file = "WGCNA_RESULT.RData",net)
load("WGCNA_RESULT.RData")

# Setting the soft power
sft=18

# Different combinations of minimum module size and deep split
for (mm in c(20,30,40,50,70,100, 150, 200, 250, 300)){
  print(mm)
  for (ds in seq(0,4)){
    print(ds)
    ds = ds; minModSize = mm; dthresh = 0.1; pam = FALSE
    networks=list()
    networks$datExpr=Proteomics_input_data
    networks$tree = hclust(1-TOM, method="average")
    networks$cut = cutreeHybrid(dendro = networks$tree, pamStage=pam, minClusterSize= minModSize, 
                                cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-TOM))
    networks$merged= mergeCloseModules(exprData= t(networks$datExpr), colors = networks$cut$labels, cutHeight=dthresh)
    networks$MEs = moduleEigengenes(t(networks$datExpr), colors=networks$merged$colors, softPower=sft)
    networks$kMEtable = signedKME(t(networks$datExpr), datME = networks$MEs$eigengenes,corFnc = "bicor")
    tryCatch(
      {
        pdf(paste("./Prot_bicor_",mm, "_DeepSplit_",ds,"sft_power",sft,".pdf",sep = ""), height = 8, width = 8)
        plotDendroAndColors(networks$tree, colors=labels2colors(networks$merged$colors), dendroLabels = F)
        dev.off()
        
      },
      error = function(e) {
        message("Too many modules to plot, skipping...")
      }
    )
    
    save(file=paste("./Prot_bicor",mm, "_DeepSplit_",ds,"sft_power",sft,".RData",sep = ""), networks)
  }
}

# Getting rid of those pooled samples
metadata_input<-metadata_input[-c(9:10),]

# module-trait plots
pdf("WGCNA_plots_module_trait_counts_norm_sft_18_bicor.pdf")
for (mm in c(20,30,40,50,70,100,150,200, 250, 300)){
  print(mm)
  for (ds in seq(0,4)){ 
    print(ds)
    load(paste("./Prot_bicor",mm, "_DeepSplit_",ds,"sft_power",sft,".RData",sep = ""))
    datExpr=networks$datExpr
    geneTree=networks$tree
    merged=networks$merged
    modules=merged$colors
    print(length(unique(modules))-1)
    MEs=networks$MEs
    kMEtable=networks$kMEtable
    table(rownames(metadata_input) == colnames(Proteomics_input_data))
    metadata_input$SampleGroup=factor(metadata_input$SampleGroup,levels = c("WT","Cul3_HET"))
    
    modTrait=data.frame()
    # Finding module trait association for each module except for the grey module
    for(i in 2:length(unique(modules))) {
      me = MEs$eigengenes[,i]
      moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
      moduleColor = labels2colors(moduleNumber)
      s = summary(lm(me ~ SampleGroup + Gender,data=metadata_input))$coefficients
      for(grp in c("Cul3_HET")){
        
        rowID = paste0("SampleGroup", grp)
        
        modTrait = rbind.fill(modTrait,
                              data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                                         beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))}
      
      modTrait$fdr=p.adjust(modTrait$p,method = "fdr")
      modTrait$signedLog10fdr = -log10(modTrait$fdr) * sign(modTrait$beta)
      modTrait$signedLog10fdr[modTrait$fdr > .05] = 0
      modTrait$text = signif(modTrait$beta, 1)
      
    }
    
    # Plotting the module-trait association
    print(ggplot(modTrait, aes(x=Module,y=Group, label=text)) +
            geom_tile(aes(fill=signedLog10fdr),color="grey60") +
            scale_fill_gradient2(low = "blue", high = "red","[beta]\nsigned\n-log10FDR\n") +
            geom_text(size=3, color="black")  +
            ggtitle(paste("soft_power_", sft,"_module_trait_module_size_", mm, "_DeepSplit_",ds, sep="")))
    
  }}
dev.off()

# GO Enrichment using gProfileR
library(gProfileR)
for (mm in c(20,30, 40, 50, 70, 100, 150, 200, 250, 300)){
  print(mm)
  for (ds in seq(0,4)){
    print(ds)
    load(paste("./Prot_bicor",mm, "_DeepSplit_",ds,"sft_power",sft,".RData",sep = ""))
    datExpr=networks$datExpr
    genes=rownames(datExpr)
    merged=networks$merged
    modules=merged$colors
    MEs = networks$MEs
    kMEtable=networks$kMEtable
    
    unique_module_labels <- sort(unique(modules))
    unique_module_labels <- unique_module_labels[unique_module_labels != 0]
    unique_module_colours <- WGCNA::labels2colors(unique_module_labels)
    
    GO_enrich<-list()
    for (i in seq_along(unique_module_labels)) {
      for (filter in c("strong", "none", "moderate")){
        moduleNumber = unique_module_labels[[i]]
        moduleColor = unique_module_colours[[i]]
        moduleGenes = genes[modules==moduleNumber]
        moduleGenes=moduleGenes[order(kMEtable[moduleGenes,paste0("kME",moduleNumber)], decreasing = T)]
        annotation<-readxl::read_xlsx("Normalized_Data.xlsx", sheet = 1)
        annotation<-as.data.frame(annotation)
        moduleGenes=annotation[match(moduleGenes, annotation$Protein_ID),]$Gene
        go = gprofiler(query=moduleGenes, 
                       organism = "mmusculus",
                       correction_method = "fdr",
                       exclude_iea = T,
                       hier_filtering = filter,
                       custom_bg = genes, 
                       src_filter = c("GO:BP","GO:MF"),
                       ordered_query = F)%>%
          arrange(p.value)%>%
          filter(overlap.size > 1)
        module_name <- paste(moduleNumber, moduleColor)
        GO_enrich[[module_name]] <- go
      }
      writexl::write_xlsx(
        GO_enrich, 
        path = paste("GOEnrich_sft18_signed_Prot_bicor_",filter,"_",ds,"_MM_",mm,".xlsx", sep = "")
        
      )
      
    }
  }
}
