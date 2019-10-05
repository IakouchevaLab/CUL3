library(tidyverse)
library(ggdendro)

load("data/source/ModuleDendrogramSource/CX_norm_counts_bicor100_DeepSplit_0sft_power18.RData")
expr <- readxl::read_xlsx("data/source/ModuleDendrogramSource/E17.5_CX_combat_norm_data_log.xlsx")

mes <- networks$MEs$eigengenes[, colnames(networks$MEs$eigengenes) != "ME0"]
me_colors <- setNames(
    WGCNA::labels2colors(as.numeric(substr(colnames(mes), 3, nchar(mes)))),
    colnames(mes)
)

me_dend <- hclust(dist(t(mes)), method = "average") %>%
    as.dendrogram() %>%
    dendro_data()
seg_coord <- segment(me_dend)
leaf_coord <- me_dend$labels %>%
    mutate(label = as.character(label)) %>%
    mutate(colour = me_colors[label])

MEDendrogram <- ggplot() +
    geom_tile(
        data = leaf_coord, mapping = aes(x = x, y = y), fill = NA, colour = NA
    ) +
    geom_segment(
        data = seg_coord,
        mapping = aes(x = x, y = y, xend = xend, yend = yend)
    ) +
    geom_point(
        data = leaf_coord,
        mapping = aes(x = x, y = y, fill = colour),
        shape = 21, size = 5, show.legend = FALSE
    ) +
    geom_text(
        data = leaf_coord,
        mapping = aes(x = x, y = y, label = label),
        hjust = 0, angle = 270, nudge_y = -0.1, size = 8
    ) +
    scale_fill_manual(values = setNames(nm = leaf_coord$colour)) +
    scale_x_continuous(expand = c(0, 0)) +
    ggnetwork::theme_blank() +
    theme(
        plot.margin = margin()
    )

me_df <- mes %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(Sample = rownames(.)) %>%
    left_join(
        readxl::read_xlsx(
            "data/source/ModuleDendrogramSource/meta_group_wgcna.xlsx",
            sheet = 1
        ), by = "Sample"
    ) %>%
    mutate(SampleGroup = factor(SampleGroup, levels = c("WT", "Cul3_HET")))
me_association <- lapply(
    colnames(me_df)[grepl("ME", colnames(me_df))],
    function(me) {
        lm(
            as.formula(paste(
                me, "~ SampleGroup + Gender"
            )),
            data = me_df
        ) %>%
            summary() %>%
            broom::tidy() %>%
            mutate(label = me) %>%
            slice(2)
    }
) %>%
    bind_rows() %>%
    mutate(FDR = p.adjust(p.value, method = "BH")) %>%
    mutate(Star = ifelse(FDR <= 0.05, "*", "")) %>%
    left_join(
        leaf_coord, by = "label"
    )

ModuleAssociation <- ggplot(
    data = me_association,
    mapping = aes(
        x = x, y = y, fill = estimate, label = Star
    )
) +
    geom_tile(colour = "black") +
    geom_text(size = 20, vjust = 0.75) +
    labs(
        y = "HET"
    ) +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    guides(
        fill = guide_colourbar(title = "Linear\nRegression\nBeta")
    ) +
    theme(
        text = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = 0, unit = "lines")
    )

sheets <- readxl::excel_sheets(
    "data/source/ASDRelevantGeneListsFromLiterature.xlsx"
)
lists <- lapply(
    setNames(
        nm = sheets[1:(length(sheets) - 1)]
    ),
    function(sheet) {
        readxl::read_xlsx(
            "data/source/ASDRelevantGeneListsFromLiterature.xlsx",
            sheet = sheet
        )[[1]]
    }
)[c(
    "CHD8Targets", "FMRPBindingTarget", "SFARISyndromicRisk12", 
    "pLI099", "SatterstromASD", "SynaptomeDBPre", "SynaptomeDBPost"
)] %>%
    setNames(
        nm = c(
            "CHD8Targets", "FMRPBindingTarget", "SFARI_S12", 
            "pLI099", "SatterstromASD", "SynaptomeDBPre", "SynaptomeDBPost"
        )
    )
module_assignments <- data.frame(
    Protein = rownames(networks$datExpr),
    Gene = toupper(expr$Gene[match(rownames(networks$datExpr), expr$Protein_ID)]),
    module_label = networks$merged$colors,
    label = paste0("ME", networks$merged$colors)
)
fisher_tests <- expand.grid(
    List = names(lists),
    module_label = sort(unique(module_assignments$module_label))
) %>%
    mutate(label = paste0("ME", module_label)) %>% 
    bind_cols(
        mapply(
            function(l, m) {
                list_genes <- toupper(lists[[l]])
                module_genes <- module_assignments %>%
                    filter(label == m) %>%
                    pull(Gene)
                all_genes <- toupper(
                    c(pull(module_assignments, Gene), unique(unlist(lists)))
                )
                a <- length(intersect(
                    module_genes, list_genes
                ))
                b <- length(intersect(
                    module_genes, all_genes[!all_genes %in% list_genes]
                ))
                c <- length(intersect(
                    list_genes, all_genes[!all_genes %in% module_genes]
                ))
                d <- length(all_genes[!all_genes %in% c(list_genes, module_genes)])
                fisher.test(
                    matrix(c(a, b, c, d), ncol = 2, nrow = 2),
                    alternative = "greater"
                ) %>%
                    broom::tidy() %>%
                    mutate(Overlap = a)
            }, .$List, .$label, SIMPLIFY = FALSE
        ) %>%
            bind_rows()
    ) %>%
    filter(module_label != 0) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    mutate(
        FDR = ifelse(estimate < 1, 1, adj.P.Val)
    ) %>%
    mutate(Star = ifelse(FDR <= 0.1, "*", "")) %>%
    mutate(
        tile_label = mapply(
            function(OR, OV, STAR) {
                paste(
                    paste0(round(OR, 2), STAR),
                    paste0("(", OV, ")"),
                    sep = "\n"
                )
            }, 
            estimate, Overlap, Star
        )
    ) %>%
    distinct() %>%
    left_join(
        leaf_coord, by = "label"
    )

ModuleListEnrichment <- ggplot(
    data = fisher_tests,
    mapping = aes(
        x = x, y = List, fill = FDR, label = tile_label
    )
) +
    geom_tile(colour = "grey") +
    geom_text(size = 8) +
    scale_fill_gradientn(
        colours = c("white", "white", "red"),
        values = c(1, 0.1, 0),
        breaks = c(1, 0.1)
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
        text = element_text(size = 25),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(t = 0, unit = "lines")
    )

ggsave(
    filename = "data/figures/E175_CX_ModuleDendroHead.pdf",
    plot = cowplot::plot_grid(
        cowplot::plot_grid(
            MEDendrogram, 
            ModuleAssociation +
                theme(
                    # plot.margin = margin(t = -25, unit = "lines"),
                    legend.position = "none"
                ), 
            ModuleListEnrichment +
                theme(
                    # plot.margin = margin(t = -40, unit = "lines"),
                    legend.position = "none"
                ),
            align = "v", axis = "lr", ncol = 1, 
            rel_heights = c(0.4, 0.075, 0.9)
        ),
        cowplot::plot_grid(
            cowplot::get_legend(
                ModuleAssociation +
                    theme(
                        # plot.margin = margin(t = -100, unit = "lines"),
                        legend.key.size = unit(3, "lines"),
                        legend.position = "bottom"
                    )
            ), 
            cowplot::get_legend(
                ModuleListEnrichment +
                    theme(
                        # plot.margin = margin(t = -40, unit = "lines"),
                        legend.position = "bottom"
                    )
            ),
            ncol = 2
        ),
        nrow = 2, rel_heights = c(1, 0.1)
    ),
    device = "pdf", width = 16, height = 16
)

