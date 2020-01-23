setwd("/Users/pramod/Desktop/Clu3/WGCNA_Cul3/E17.5_WGCNA_CX/WGCNA_Proteomics/KME_table/")
library(tidyverse)
library(UniProt.ws)
library(WGCNA)
library(ggrepel)
enableWGCNAThreads()
allowWGCNAThreads()


################################################################################
# Prepare data                                                                 #
################################################################################

load("./CX_E17.5_norm_counts_bicor100_DeepSplit_0sft_power18.RData")
net_E17.5_CX<- networks
assign_E17.5_CX <- readxl::read_xlsx(
    "E17.5_CX_WGCNA_Proteomics_KME_Table.xlsx",
  sheet = 1
)

load("CX_P7_norm_counts_bicor100_DeepSplit_0sft_power18.RData")
net_P7_CX <- networks
assign_P7_CX <- readxl::read_xlsx(
  "P7_CX_WGCNA_Proteomics_KME_Table.xlsx",
  sheet = 1
)

load("CX_P35_norm_counts_bicor20_DeepSplit_0sft_power18.RData")
net_P35_CX <- networks
assign_P35_CX <- readxl::read_xlsx(
  "P35_CX_WGCNA_Proteomics_KME_Table.xlsx",
  sheet = 1
)



load("./CB_E17.5_norm_counts_bicor40_DeepSplit_0sft_power22.RData")
net_E17.5_CB<- networks
assign_E17.5_CB <- readxl::read_xlsx(
  "E17.5_CB_WGCNA_Proteomics_KME_Table.xlsx",
  sheet = 1
)

load("CB_P7_norm_counts_bicor150_DeepSplit_0sft_power26.RData")
net_P7_CB <- networks
assign_P7_CB <- readxl::read_xlsx(
  "P7_CB_WGCNA_Proteomics_KME_Table.xlsx",
  sheet = 1
)

load("CB_P35_norm_counts_bicor40_DeepSplit_2sft_power12.RData")
net_P35_CB <- networks
assign_P35_CB <- readxl::read_xlsx(
  "P35_CB_WGCNA_Proteomics_KME_Table.xlsx",
  sheet = 1)
################################################################################
# Module Preservation                                                          #
################################################################################


protein_intersection <-  Reduce(intersect, list(rownames(net_E17.5_CX$datExpr), rownames(net_P7_CX$datExpr), 
                                                rownames(net_P35_CX$datExpr), rownames(net_E17.5_CB$datExpr), 
                                                rownames(net_P7_CB$datExpr), rownames(net_P35_CB$datExpr)))



multiexpr <- list(
  "E17.5_CX" = list(data = t(net_E17.5_CX$datExpr[protein_intersection, ])),
  "P7_CX" = list(data = t(net_P7_CX$datExpr[protein_intersection, ])),
  "P35_CX" = list(data = t(net_P35_CX$datExpr[protein_intersection, ])),
  "E17.5_CB" = list(data = t(net_E17.5_CB$datExpr[protein_intersection, ])),
  "P7_CB" = list(data = t(net_P7_CB$datExpr[protein_intersection, ])),
  "P35_CB" = list(data = t(net_P35_CB$datExpr[protein_intersection, ]))
  
)
checkSets(multiexpr)

multicolor <- list(
  "E17.5_CX" = assign_E17.5_CX$Module[match(protein_intersection, assign_E17.5_CX$Accessions)],
  "P7_CX" = assign_P7_CX$Module[match(protein_intersection, assign_P7_CX$Accessions)],
  "P35_CX" = assign_P35_CX$Module[match(protein_intersection, assign_P35_CX$Accessions)],
  "E17.5_CB" = assign_E17.5_CB$Module[match(protein_intersection, assign_E17.5_CB$Accessions)],
  "P7_CB" = assign_P7_CB$Module[match(protein_intersection, assign_P7_CB$Accessions)],
  "P35_CB" = assign_P35_CB$Module[match(protein_intersection, assign_P35_CB$Accessions)]
)

mp <- modulePreservation(
  multiData = multiexpr,
  multiColor = multicolor,
  referenceNetworks = c(1)
)


mp_z <- bind_rows(
  mp$preservation$Z$ref.E17.5_CX$inColumnsAlsoPresentIn.P7_CX %>%
    mutate(module_colour = rownames(.)) %>%
    mutate(Reference = 1, Test = 2) %>%
    mutate(
      strip_label = "Reference: E17.5_CX; Test: P7_CX"
    ),
  mp$preservation$Z$ref.E17.5_CX$inColumnsAlsoPresentIn.P35_CX %>%
    mutate(module_colour = rownames(.)) %>%
    mutate(Reference = 1, Test = 2) %>%
    mutate(
      strip_label = "Reference: E17.5_CX; Test: P35_CX"
    ),
  mp$preservation$Z$ref.E17.5_CX$inColumnsAlsoPresentIn.E17.5_CB %>%
    mutate(module_colour = rownames(.)) %>%
    mutate(Reference = 1, Test = 2) %>%
    mutate(
      strip_label = "Reference: E17.5_CX; Test: E17.5_CB"
    ),
  mp$preservation$Z$ref.E17.5_CX$inColumnsAlsoPresentIn.P7_CB %>%
    mutate(module_colour = rownames(.)) %>%
    mutate(Reference = 1, Test = 2) %>%
    mutate(
      strip_label = "Reference: E17.5_CX; Test: P7_CB"
    ),
  mp$preservation$Z$ref.E17.5_CX$inColumnsAlsoPresentIn.P35_CB %>%
    mutate(module_colour = rownames(.)) %>%
    mutate(Reference = 1, Test = 2) %>%
    mutate(
      strip_label = "Reference: E17.5_CX; Test: P35_CB"
    )
)

mp_plot <- ggplot(
  data = mp_z_A %>%
    filter(module_colour != c("gold")),
  mapping = aes(
    x = moduleSize, y = Zsummary.pres,
    fill = module_colour, label = module_colour
  )
) +
  facet_wrap(~ strip_label, scales = "free") +
  geom_point(shape = 21, colour = "black", size = 20, show.legend = FALSE) +
  geom_hline(yintercept = 2, colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 10, colour = "green", linetype = "dashed") +
  geom_text_repel(
    data = mp_z_A %>%
      filter(module_colour != c("gold")),
    size = 15
  )  +
  labs(
    x = "Module Size", y = "Preservation Z-Summary"
  ) +
  scale_fill_manual(
    values = setNames(nm = mp_z_A$module_colour)
  ) +
  scale_x_continuous(trans = "log1p", breaks = c(10, 50, 100, 500, 1000)) +
  theme_bw() +
  theme(
    text = element_text(size = 50)
  )

ggsave(
  filename = "./ModulePreservation_CX_overlap_all_other.pdf",
  plot = mp_plot,
  device = "pdf", width = 70, height = 40,
  limitsize = FALSE,
  useDingbats = FALSE
)


