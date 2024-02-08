library(tidyverse)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

gsea_results <- read_rds("data_generated/gsea.rds")

selected_comparisons <- c(
  "Resting" = "Resting_Target_vs_NTC",
  "Tgfb1" = "Tgfb1_Target_vs_NTC",
  "Il1b" = "Il1b_Target_vs_NTC"
)



# Plot data ---------------------------------------------------------------

## Check NES ----

ggplot(gsea_results, aes(NES, -log10(padj))) +
  geom_point(alpha = 0.25)


## Plot selected terms ----

selected_terms <- tribble(
  ~db, ~pathway, ~term_display,
  "MSigDB_Hallmark_2020",
  "TNF-alpha Signaling via NF-kB",
  "TNF-alpha Signaling via NF-kB (MSigDB)",
  "MSigDB_Hallmark_2020",
  "Inflammatory Response",
  "Inflammatory Response (MSigDB)",
  "KEGG_2019_Mouse",
  "JAK-STAT signaling pathway",
  "JAK-STAT signaling pathway (KEGG)",
  "KEGG_2019_Mouse",
  "Hedgehog signaling pathway",
  "Hedgehog signaling pathway (KEGG)",
  "WikiPathways_2019_Mouse",
  "Striated Muscle Contraction WP216",
  "Striated Muscle Contraction (WikiPathways)",
  "MSigDB_Hallmark_2020",
  "Oxidative Phosphorylation",
  "Oxidative Phosphorylation (MSigDB)",
  "Reactome_2022",
  "Mitochondrial Translation R-HSA-5368287",
  "Mitochondrial Translation (Reactome)",
  "WikiPathways_2019_Mouse",
  "Retinol metabolism WP1259",
  "Retinol metabolism (WikiPathways)",
  "MSigDB_Hallmark_2020",
  "G2-M Checkpoint",
  "G2-M Checkpoint (MSigDB)"
)

selected_targets <- c(
  "Tgfbr1", "Smad2", "Smad3", "Chd4", "Kat5",
  "Dmap1", "Srcap", "Smarca4", "Brd9", "Brd7",
  "Arid2", "Pbrm1", "Kat8", "Hcfc1", "Wdr82",
  "Paxip1", "Kmt2a", "Egr2", "Rest", "Rnf40"
)

plot_terms <- function(terms, targets) {
  terms <- 
    terms %>% 
    mutate(term_display = as_factor(term_display))
  
  # prepare data
  data_vis <- 
    gsea_results %>% 
    filter(
      comparison %in% selected_comparisons,
      group %in% targets
    ) %>% 
    left_join(terms, by = join_by(db, pathway)) %>% 
    filter(!is.na(term_display)) %>% 
    mutate(
      group = fct_relevel(group, !!!targets),
      comparison = 
        comparison %>% 
        fct_recode(!!!selected_comparisons) %>% 
        fct_relevel(!!!names(selected_comparisons))
    )
  
  # set limits for colorbar
  color_limit <- max(abs(data_vis$NES), na.rm = TRUE)
  
  # make plot  
  ggplot(data_vis, aes(comparison, term_display, size = -log10(padj))) +
    geom_point(aes(color = NES)) +
    xlab("Condition") +
    ylab("Term (database)") +
    scale_color_distiller(
      name = "Normalized\nenrichment\nscore",
      palette = "RdBu",
      direction = -1,
      limits = c(-color_limit, color_limit),
      breaks = c(-color_limit, 0, color_limit),
      labels = function(x) round(x, 2),
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_size_area(
      name = "-log10 padj",
      max_size = 3,
      limits = c(0, 10),
      breaks = c(0, 5, 10),
      labels = c("0", "5", "10 or higher"),
      oob = scales::oob_squish
    )  +
    coord_fixed() +
    ggtitle("Target vs NTC") +
    facet_wrap(vars(group), nrow = 1) +
    theme_pub(TRUE) +
    theme(
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.margin = margin(),
      panel.grid = element_blank(),
      panel.spacing = unit(-.5, "pt")
    )
}

plot_terms(selected_terms, selected_targets)
ggsave_default("4d_gsea", type = "pdf", width = 210)


## Plot GSEA per database ----

plot_gsea <- function(db = "MSigDB_Hallmark_2020",
                      comparisons = c(
                        "Il1b_Target_vs_NTC",   
                        "Resting_Target_vs_NTC",
                        "Tgfb1_Target_vs_NTC"
                      ),
                      groups = NULL,
                      top_n_positive = 10L,
                      top_n_negative = 10L,
                      terms = NULL,
                      max_p_adj = 0.05,
                      min_abs_NES = 1,
                      source_data_filename = NULL,
                      source_data_sheet = NULL) {
  # prepare data
  data <-
    gsea_results %>%
    filter(db == {{db}})
  
  if (!is.null(comparisons))
    data <- 
      data %>% 
      filter(comparison %in% {{comparisons}})
  
  if (!is.null(groups))
    data <- 
      data %>% 
      filter(group %in% {{groups}})
  
  if (is.null(terms)) {
    data_top_terms <-
      data %>% 
      filter(
        padj <= max_p_adj,
        abs(NES) >= min_abs_NES
      ) %>% 
      group_by(comparison, group)
    
    top_terms_pos <- 
      data_top_terms %>% 
      slice_max(n = top_n_positive, order_by = NES, with_ties = FALSE) %>%
      pull(pathway) %>%
      unique()
    
    top_terms_neg <- 
      data_top_terms %>% 
      slice_min(n = top_n_negative, order_by = NES, with_ties = FALSE) %>%
      pull(pathway) %>%
      unique()
    
    terms <- c(top_terms_pos, top_terms_neg)
    # info(terms)
  }
  
  data_vis <- 
    data %>% 
    filter(pathway %in% terms) %>% 
    mutate(
      pathway =
        as_factor(pathway) %>%
        fct_reorder(NES * -log10(padj), sum, na.rm = TRUE),
      # group = factor(group),
      comparison = 
        comparison %>% 
        fct_recode(!!!selected_comparisons) %>% 
        fct_relevel(!!!names(selected_comparisons))
    )
  
  if (nlevels(data_vis$pathway) > 5) {
    horizontal_grid <-
      geom_hline(
        yintercept = seq(6, nlevels(data_vis$pathway), 5) - 0.5,
        linewidth = BASE_LINEWIDTH,
        color = "grey92"
      )
  } else {
    horizontal_grid <- NULL  
  }
  
  color_limit <- max(abs(gsea_results$NES), na.rm = TRUE)
  
  # make plot  
  ggplot(data_vis, aes(comparison, pathway, size = -log10(padj))) +
    scale_y_discrete() +
    horizontal_grid +    
    geom_point(aes(color = NES)) +
    xlab("condition") +
    ylab(str_glue("{str_replace_all(db, '_', ' ')} gene set")) +
    scale_color_distiller(
      name = "normalized\nenrichment\nscore",
      palette = "RdBu",
      direction = -1,
      limits = c(-color_limit, color_limit),
      breaks = c(-color_limit, 0, color_limit),
      labels = function(x) round(x, 2),
      guide = guide_colorbar(
        barheight = unit(15, "mm"),
        barwidth = unit(2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_size_area(
      name = "-log10 padj",
      max_size = 3,
      limits = c(0, 10),
      breaks = c(0, 5, 10),
      labels = c("0", "5", "10 or higher"),
      oob = scales::oob_squish
    )  +
    coord_fixed() +
    ggtitle("target vs NTC") +
    facet_wrap(vars(group), nrow = 1) +
    theme_pub(TRUE) +
    theme(
      legend.key.height = unit(3, "mm"),
      legend.key.width = unit(3, "mm"),
      legend.margin = margin(),
      panel.grid = element_blank(),
      panel.spacing = unit(-.5, "pt")
    )
}

# plot_gsea("KEGG_2019_Mouse")
# ggsave_default("gsea_KEGG_2019_Mouse", type = "pdf")
# 
# plot_gsea("MSigDB_Hallmark_2020")
# ggsave_default("gsea_MSigDB_Hallmark_2020", type = "pdf")
# 
# plot_gsea("NCI-Nature_2016")
# ggsave_default("gsea_NCI-Nature_2016", type = "pdf")
# 
# plot_gsea("Reactome_2022")
# ggsave_default("gsea_Reactome_2022", type = "pdf")
# 
# plot_gsea("WikiPathways_2019_Mouse")
# ggsave_default("gsea_WikiPathways_2019_Mouse", type = "pdf")
# 
# plot_gsea("ChEA_2022", top_n_positive = 5, top_n_negative = 5)
# ggsave_default("gsea_ChEA_2022", type = "pdf")
# 
# plot_gsea("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
# ggsave_default("gsea_ENCODE_and_ChEA_Consensus", type = "pdf")
# 
# plot_gsea("ENCODE_TF_ChIP-seq_2015", top_n_positive = 5, top_n_negative = 5)
# ggsave_default("gsea_ENCODE_TF_ChIP-seq_2015", type = "pdf")
# 
# plot_gsea("TRANSFAC_and_JASPAR_PWMs")
# ggsave_default("gsea_TRANSFAC_and_JASPAR_PWMs", type = "pdf")
# 
# plot_gsea("TRRUST_Transcription_Factors_2019")
# ggsave_default("gsea_TRRUST_Transcription_Factors_2019", type = "pdf")
