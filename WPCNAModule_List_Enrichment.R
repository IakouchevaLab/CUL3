library(tidyverse)
library(ggdendro)

# Data load
gene_lists <- lapply(
    setNames(nm = readxl::excel_sheets('src/ASDRelevantGeneListsFromLiterature.xlsx')[-13]),
    function(sheet) {
        readxl::read_xlsx(
            'src/ASDRelevantGeneListsFromLiterature.xlsx',
            sheet = sheet
        )[, 1] %>%
            distinct()
    }
)
gene_lists <- gene_lists[!names(gene_lists) %in% c(
    "EichlerDNM_LGD_MIS_ASD_BD_SCZ", 
    "EichlerDNM_LGD_ASD_BD_SCZ",
    "VulnerableASD",
    "ASDSanders65"
)]

eigengene_files <- list.files(
    'src', pattern = "*eigengene*", full.names = TRUE
)
names(eigengene_files) <- lapply(eigengene_files, function(fn) {
    paste(
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 1],
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 2],
        sep = "_"
    )
})
eigengene_data <- lapply(
    eigengene_files,
    function(fn) {
        readxl::read_xlsx(fn, sheet = 1) %>%
            mutate(filename = basename(fn)) %>%
            mutate(
                Reg = str_split(filename, pattern = "_", simplify = TRUE)[, 1],
                Age = str_split(filename, pattern = "_", simplify = TRUE)[, 2]
            ) %>%
            mutate(Sample = str_split(Samplename, "_", simplify = TRUE)[, 1])
    }
)
module_membership_files <- list.files(
    'src', pattern = "*module_genes*", full.names = TRUE
)
names(module_membership_files) <- lapply(module_membership_files, function(fn) {
    paste(
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 4],
        str_split(basename(fn), pattern = "_", simplify = TRUE)[, 3],
        sep = "_"
    )
})
module_membership_data <- lapply(
    module_membership_files,
    function(fn) {
        lapply(
            setNames(nm = readxl::excel_sheets(fn)),
            function(sheet) {
                readxl::read_xlsx(fn, sheet = sheet) %>%
                    mutate(filename = basename(fn)) %>%
                    mutate(
                        Reg = str_split(filename, pattern = "_", simplify = TRUE)[, 4],
                        Age = str_split(filename, pattern = "_", simplify = TRUE)[, 3]
                    ) %>%
                    mutate(label = sheet)
            }
        ) %>%
            bind_rows()
    }
)
metadata <- lapply(
    Sys.glob("src/meta_group*\\.xlsx"),
    function(fn) {
        readxl::read_xlsx(fn, sheet = 1) %>%
            mutate(filename = fn) %>%
            mutate(
                Reg = str_split(filename, pattern = "_", simplify = TRUE)[, 4] %>%
                    str_replace(".xlsx", ""),
                Age = str_split(filename, pattern = "_", simplify = TRUE)[, 3]
            )
    }
) %>%
    bind_rows() %>%
    mutate(Genotype = sapply(SampleGroup, function(x) {
        if (x == "Cul3 HET") { return("HET") } else { return(x) }
    })) %>%
    filter(Genotype %in% c("WT", "HET")) %>%
    mutate(Genotype = relevel(factor(Genotype), ref = "WT"))

# Define Region-Age combinations
reg_age <- metadata[, c("Reg", "Age")] %>%
    arrange(Reg, Age) %>%
    distinct() %>%
    mutate(RegAge = paste(Reg, Age, sep = "_"))

all_plots <- list()
for (i in seq(1, nrow(reg_age))) {
    this_reg <- as.character(reg_age[i, "Reg"])
    this_age <- as.character(reg_age[i, "Age"])
    this_reg_age <- as.character(reg_age[i, "RegAge"])
    this_eigengene <- eigengene_data[[this_reg_age]] %>%
        inner_join(
            metadata, 
            by = c("Sample" = "SampleName", "Reg", "Age")
        )
    this_module_membership <- module_membership_data[[this_reg_age]]
    all_genes <- unique(c(toupper(this_module_membership), toupper(unlist(gene_lists))))
    # Calculate all data
    # Module-trait association
    module_trait <- lapply(
        setNames(nm = colnames(dplyr::select(this_eigengene, starts_with("ME")))),
        function(module) {
            this_formula <- paste(module, "~ Genotype")
            lm(as.formula(this_formula), data = this_eigengene) %>%
                broom::tidy() %>%
                mutate(label = module) %>%
                mutate(fdr = p.adjust(p.value, method = "BH")) %>%
                mutate(Genotype = c("WT", "HET")) %>%
                filter(Genotype == "HET")
        }
    ) %>%
        bind_rows() %>%
        mutate(Star = ifelse(fdr <= 0.05, "*", ""))
    # List enrichment
    list_enrichment <- lapply(
        setNames(nm = unique(this_module_membership$label)),
        function(ml) {
            lapply(
                names(gene_lists),
                function(gl) {
                    this_list <- toupper(pull(gene_lists[[gl]], gene_symbol))
                    module_list <- toupper(pull(filter(this_module_membership, label == ml), Gene_name))
                    contingency <- matrix(
                        c(
                            length(intersect(this_list, module_list)),
                            length(intersect(this_list, all_genes[!all_genes %in% module_list])),
                            length(intersect(all_genes[!all_genes %in% this_list], module_list)),
                            length(all_genes[!all_genes %in% c(this_list, module_list)])
                        ),
                        ncol = 2, nrow = 2, byrow = TRUE
                    )
                    fisher.test(contingency, alternative = "greater") %>%
                        broom::tidy() %>%
                        mutate(gene_list = gl) %>%
                        mutate(label = str_replace(ml, "M", "ME"))
                }
            ) %>%
                bind_rows()
        }
    ) %>%
        bind_rows() %>%
        mutate(fdr = p.adjust(p.value, method = "BH")) %>%
        mutate(Star = ifelse(fdr <= 0.05, "*", ""))
    selected_modules <- unique(c(
        pull(filter(module_trait, Star == "*"), label),
        pull(filter(list_enrichment, Star == "*"), label)
    ))
    selected_modules <- selected_modules[selected_modules != "ME0"]
    # Dendrogram
    cluster_MEs <- this_eigengene[, colnames(this_eigengene) %in% selected_modules] %>%
        dplyr::select(starts_with("ME")) %>%
        t() %>%
        as.matrix() %>%
        dist() %>%
        hclust()
    dendrogram_MEs <- as.dendrogram(cluster_MEs)
    dendrogram_data_MEs <- dendro_data(dendrogram_MEs)
    dendrogram_segments_MEs <- segment(dendrogram_data_MEs)
    # This data will be used to align ALL plots
    dendrogram_labels_MEs <- dendrogram_data_MEs$labels %>%
        mutate(module_label = as.numeric(str_replace(label, "ME", ""))) %>%
        mutate(module_colour = WGCNA::labels2colors(module_label))
    dendrogram_plt <- ggplot() +
        coord_equal() +
        geom_segment(
            data = dendrogram_segments_MEs,
            mapping = aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_point(
            data = dendrogram_labels_MEs,
            mapping = aes(x = x, y = y, fill = module_colour),
            shape = 21, show.legend = FALSE, size = 5
        ) +
        geom_text(
            data = dendrogram_labels_MEs,
            mapping = aes(x = x, y = y, label = label), 
            size = 5, nudge_y = -0.5, angle = -90, hjust = 0
        ) +
        geom_tile(
            data = dendrogram_labels_MEs,
            mapping = aes(x = x, y = y),
            alpha = 0
        ) +
        scale_fill_manual(values = setNames(nm = dendrogram_labels_MEs$module_colour)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.75, 0.001))) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    dendrogram_plt
    module_trait_plot_data <- left_join(
        module_trait, dendrogram_labels_MEs, by = "label"
    ) %>%
        filter(label != "ME0")
    module_trait_plt <- ggplot() +
        coord_equal() +
        geom_tile(
            data = module_trait_plot_data,
            mapping = aes(x = x, y = Genotype, fill = estimate),
            colour = "grey"
        ) +
        geom_text(
            data = module_trait_plot_data,
            mapping = aes(x = x, y = Genotype, label = Star),
            size = 10, vjust = 0.75
        ) +
        scale_fill_gradient2(
            low = "lightblue", mid = "white", high = "firebrick1", 
            midpoint = 0, limits = c(-1, 1)
        ) +
        guides(
            fill = guide_colourbar(title = "Linear Regression Beta")
        ) +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    module_trait_plt
    list_enrichment_plot_data <- left_join(
        list_enrichment, dendrogram_labels_MEs, by = "label"
    ) %>%
        filter(label != "ME0")
    list_enrichment_plt <- ggplot() +
        coord_equal() +
        geom_tile(
            data = list_enrichment_plot_data %>%
                filter(Star == "*"),
            mapping = aes(x = x, y = gene_list, fill = -log10(fdr)),
            colour = "grey"
        ) +
        geom_tile(
            data = list_enrichment_plot_data %>%
                filter(Star == "") %>%
                filter(label %in% selected_modules),
            mapping = aes(x = x, y = gene_list),
            colour = "grey", fill = "white"
        ) +
        # geom_tile(
        #     data = list_enrichment_plot_data,
        #     mapping = aes(x = x, y = gene_list, fill = -log10(fdr)),
        #     colour = "grey"
        # ) +
        geom_text(
            data = list_enrichment_plot_data,
            mapping = aes(x = x, y = gene_list, label = Star),
            size = 10, vjust = 0.75
        ) +
        scale_fill_gradient(high = "darkorchid", low = "white") +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank()
        )
    list_enrichment_plt
    writexl::write_xlsx(
        x = list("ModuleTraitAssociation" = module_trait_plot_data, "ListEnrichment" = list_enrichment_plot_data),
        path = paste0("data/WPCNA_TraitAssoc_ListEnr_", this_reg_age, ".xlsx")
    )
    ggsave(
        filename = paste0("data/figures/", this_reg_age, "_Purple_NegLog10FDR_WPCNA.pdf"),
        plot = egg::ggarrange(
            dendrogram_plt,
            module_trait_plt,
            list_enrichment_plt,
            ncol = 1
        ),
        width = 12, height = 12, device = "pdf"
    )
}
