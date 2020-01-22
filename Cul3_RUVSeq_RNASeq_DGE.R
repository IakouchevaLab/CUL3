library(limma)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(ffpe)
library(RColorBrewer)
library(tidyverse)


# Function Name: concise_RUVSeq
# 
# Parameters: 
#         metatable_name: The file name of the metadata sheet, which contains Samples' ID, Gender,
#                         Genotype (A, B) and Brain Regions (X,Y,Z).
#         period: The period of the input sample.
#         region: The region of the input sample.
#         RUV_input_filename: count matrix after mapping to mouse genome.
#         negctrl_filename: The file name of the file which contains negative control genes.
#         num_k: The number of differnent k's that the user want to try on. (Ex. If the input is 3, this 
#               function will run RUVSeq for k = 1,2,3)
#
concise_RUVSeq <- function(metatable_name, period, region, RUV_input_filename, 
                            negctrl_filename, num_k){
  
  dir.create("./RUVSeq_Plots")
  dir.create("./RUVSeq_Output_files")
  
  metatable <- read.table(metatable_name, header = T)
  rownames(metatable) <- metatable$Sample ; metatable$Sample=NULL
  metatable_region <- metatable[metatable$Region == region,]
  
  Ruv_input=read.table(RUV_input_filename, header = T, row.names = 1)
  negControls <- read.table(negctrl_filename, sep = "\t", header = TRUE, as.is = TRUE)
  
  Region_Geno <- as.factor(rep(c("A", "B"), each=6))
  names(Region_Geno) <- colnames(Ruv_input)
  
  Region_Geno_new <- factor(Region_Geno, levels = c("B", "A"))
  
  Ruv_input <- as.matrix(Ruv_input)
  
  neg_Con <- intersect(negControls[,2], rownames(Ruv_input))

  colors <- brewer.pal(9, "Set1")
  colLib <- colors[Region_Geno]
  
  uq_region <- betweenLaneNormalization(Ruv_input, which = "upper")
  pdf(paste0("./RUVSeq_Plots/","RUVSeq_uq_",region,".pdf"), width = 50, height = 50)
  plotRLE(uq_region, col = colLib, outline=FALSE, las = 3, ylim = c(-.2, .2), 
          ylab = "Relative log Expression", cex.axis = 0.6, cex.lab = 1)
  plotPCA(uq_region, col = colLib, cex=1,cex.axis = 1.0, 
          cex.lab = 1, xlim = c(-.6, .9), ylim = c(-.7,.6) )
  dev.off()
  
  groups_region_new <- matrix(data = c(1:6, 7:12), nrow=2, byrow = TRUE)

  S_Regions <- list("uq_region" = uq_region)
  for (i in seq(num_k)) {
    S_Region <- RUVs(uq_region, neg_Con, k=i, groups_region_new)
    S_Regions <- append(S_Regions,list(S_Region))
    pdf(paste0("./RUVSeq_Plots/","RUVSeq_",period,"_",region,"_k_",i,".pdf"), width = 50, height = 50)
    plotRLE(S_Region$normalizedCounts, col = colLib,outline = FALSE,las = 3,ylim = c(-.2, .2), 
            ylab = "Relative log Expression", cex.axis = 0.6, cex.lab = 1 )
    plotPCA(S_Region$normalizedCounts, col = colLib, cex = 1,cex.axis = 1.0, 
            cex.lab = 1, xlim = c(-.6, .9), ylim = c(-.7,.6) )
    dev.off()
  
  }
  S_Regions <- set_names(S_Regions, nm = c("uq_region", seq(num_k)))
  
  options <- expand.grid(names(S_Regions),c(TRUE,FALSE))
  
  for (i in seq(length(options[[1]]))){
    
    if (options[i,][[2]] & options[i,][[1]] != "uq_region"){
      design_new <- model.matrix(~Region_Geno_new + S_Regions[[as.character(options[i,][[1]])]]$W + 
                                   metatable_region$Sex )
    } else if(options[i,][[2]] & options[i,][[1]] == "uq_region") {
      design_new <- model.matrix(~Region_Geno_new + metatable_region$Sex)
    } else if(options[i,][[1]] != "uq_region"){
      design_new <- model.matrix(~Region_Geno_new + S_Regions[[as.character(options[i,][[1]])]]$W)
    } else {
      design_new <- model.matrix(~Region_Geno_new)
    }
    
    y <- DGEList(counts=Ruv_input, group=Region_Geno_new) %>% 
      calcNormFactors(., method="upperquartile") %>% 
      estimateGLMCommonDisp(., design_new, verbose=TRUE) %>% 
      estimateGLMTagwiseDisp(., design_new)
    
    fit <- glmFit(y, design_new)
    lrt <- glmLRT(fit, coef=2)
    topFC <- topTags(lrt, n=Inf)$table
    write.table(topFC, file = if_else(options[i,][[2]],
                                        true = paste0("./RUVSeq_Output_files/",period, "_", region, 
                                                      "_", options[i,][[1]],"_S_", "edgeR_results_level.txt"),
                                       false = paste0("./RUVSeq_Output_files/",period, "_", region, 
                                                      "_", options[i,][[1]], "_edgeR_results_level.txt")), 
                sep="\t")
    
    
    pdf(paste0("./RUVSeq_Plots/","RUVSeq_",period,"_",region, options[i,][[1]],"p_val_FDR.pdf"), 
        width = 50, height = 50)
    hist(topFC$PValue, main="", xlab="p-value", breaks=100, ylim=c(0, 1400))
    hist(topFC$FDR, main="", xlab="FDR", breaks=100, ylim=c(0, 1400))
    
    plot(topFC[,1], -log10(topFC$PValue), pch=20, col="gray", cex=.5, 
         ylab="-log10(p-value)", xlab="log2(B/A)", ylim=c(0, 85), xlim=c(-2, 4), cex.lab=1, cex.axis=1)
    de <- rownames(topFC[topFC$FDR<=0.01,])
    points(topFC[de,1], -log10(topFC[de, "PValue"]), pch=20, col=colors[2], cex=1)
    dev.off()
    
  }
  
}






