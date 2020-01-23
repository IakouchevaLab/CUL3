library(tidyverse)
source("plotting.R")

# Reading the sheet names from the literature lists excel file
sheets <- readxl::excel_sheets(
  "Literature_lists.xlsx"
)
lists <- lapply(
  setNames(
    nm = sheets[1:(length(sheets) - 1)]
  ),
  function(sheet) {
    readxl::read_xlsx(
      "Literature_lists.xlsx",
      sheet = sheet
    )[[1]]
  }
)[c(
  "CHD8Targets", "FMRPBindingTarget", "SFARISyndromicRisk12", 
  "pLI099", "SatterstromASD", "SynaptomeDBPre", "SynaptomeDBPost"
)] %>%
  lapply(., unique)

# List DGE files
deg_fn <- setNames(
  list.files("./DGE", full.names = TRUE),
  list.files("./DGE/")
)

# List DPE files
dep_fn <- setNames(
  list.files("./DPE", full.names = TRUE),
  list.files("./DPE/")
)

# Extracting Ensembl_ID, Gene_name and FDR value from each
# DGE and DPE list. Add two new columns to indicate the region
# for each gene and which dataset do they come from
de_all <- bind_rows(
  lapply(
    setNames(nm = names(deg_fn)),
    function(fn) {
      period <- str_split(string = fn, pattern = "_", simplify = TRUE)[[1]]
      region <- gsub(
        ".xlsx", "",
        str_split(string = fn, pattern = "_", simplify = TRUE)[[2]]
      )
      print(fn)
      readxl::read_xlsx(
        deg_fn[[fn]], sheet = 1
      ) %>%
        dplyr::select(
          Ensembl_ID...1, Gene_name, FDR
        ) %>%
        rename(replace = c("Ensembl_ID...1" = "ID", "FDR" = "P")) %>%
        mutate(region = region)
    }
  ) %>%
    bind_rows() %>%
    mutate(Data = "DEG"),
  lapply(
    setNames(nm = names(dep_fn)),
    function(fn) {
      period <- str_split(string = fn, pattern = "_", simplify = TRUE)[[1]]
      region <- gsub(
        ".xlsx", "",
        str_split(string = fn, pattern = "_", simplify = TRUE)[[2]]
      )
      print(fn)
      readxl::read_xlsx(
        dep_fn[[fn]], sheet = 1
      ) %>%
        dplyr::select(
          Protein_ID, GeneSymbol, adj.P.Val
        ) %>%
        rename(replace = c("Protein_ID" = "ID",
                           "adj.P.Val" = "P",
                           "GeneSymbol" = "Gene_name")) %>%
        mutate(region = region)
    }
  ) %>%
    bind_rows() %>%
    mutate(Data = "DEP")
)

# Create a Data Frame to contains all possible combinations
# between target lists and literature lists which could happen 
# in the further test
fisher_tests <- de_all %>%
  distinct(region, Data) %>% {
    lapply(
      names(lists),
      function(lis) {
        mutate(., List = lis)
      }
    ) %>%
      bind_rows()
  }
fisher_tests <- lapply(
  names(lists),
  function(lis) {
    mutate(fisher_tests, List = lis)
  }
) %>%
  bind_rows()

# Performing Fisher Test and column bind the result to the 
# previously built fisher_test data frame
fisher_results <- fisher_tests %>%
  bind_cols(
    ., 
    apply(., 1, function(param) {
      this_dat <- filter(
        de_all, region == param[["region"]] &
          Data == param[["Data"]]
      )
      this_list <- param[["List"]]
      list_genes <- toupper(lists[[this_list]])
      all_genes <- toupper(c(pull(this_dat, Gene_name), unique(unlist(lists))))
      deg <- toupper(pull(filter(this_dat, P <= 0.1), Gene_name))
      a <- length(intersect(deg, list_genes))
      b <- length(intersect(
        deg, all_genes[!all_genes %in% list_genes]
      ))
      c <- length(intersect(
        list_genes, all_genes[!all_genes %in% deg]
      ))
      d <- length(all_genes[!all_genes %in% c(list_genes, deg)])
      fisher.test(
        matrix(c(a, b, c, d), ncol = 2, nrow = 2),
        alternative = "greater"
      ) %>%
        broom::tidy() %>%
        mutate(Overlap = a)
    }) %>%
      bind_rows()
  ) %>%
  filter(List != "CHD8Targets") %>%
  mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
  mutate(
    FDR = ifelse(estimate < 1, 1, adj.P.Val)
  ) %>%
  mutate(Star = ifelse(FDR <= 0.1, "*", "")) %>%
  mutate(
    region = factor(region, levels = c("CX", "CB", "HIP"))
  ) %>%
  mutate(
    label = mapply(
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
  mutate("Signed_log10_FDR" = -log10(FDR))


# Plotting
enr_all <- ggplot() +
  facet_nested(. ~ Data + region, scales = "free") +
  geom_tile(
    data = fisher_results %>%
      filter(FDR > 0.1),
    mapping = aes(x = "", y = List),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = fisher_results %>%
      filter(FDR <= 0.1),
    mapping = aes(x = "", y = List, fill = Signed_log10_FDR),
    colour = "black"
  ) +
  geom_text(
    data = fisher_results,
    mapping = aes(x = "", y = List, label = label),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


# Saving the plot
pdf("./DE_ListEnrichment_Fisher_DGE_DPE_union_by_region_FDR.pdf", width = 16, height = 9)
invisible(
  lapply(
    list(
      enr_all
    ),
    function(p) {
      p <- p + 
        theme(
          text = element_text(size = 30)
        )
      print(p)
    }
  )
)
dev.off()


