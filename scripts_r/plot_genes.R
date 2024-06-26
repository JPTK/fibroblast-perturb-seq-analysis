library(tidyverse)
library(ggrepel)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

dge <- read_rds("data_generated/dge.rds")



# Plot data ---------------------------------------------------------------

labeled_genes <- c(
  "Acta2",
  "Myl9",
  "Pdlim1",
  "Myh9",
  "Palld",
  "Tagln",
  "Tnc",
  "Ltbp2",
  "Hspg2",
  "Col1a1",
  "Col1a2",
  "Lox",
  "Ddah1",
  "Hbegf",
  "Ccn2",
  "Egr2",
  "Meox1",
  "Pdgfra",
  "Mgp",
  "Spon2",
  "Cthrc1",
  "Postn",
  "Cygb",
  "Vcan"
)



plot_volcano <- function(genes) {
  data_vis <- 
    dge %>% 
    filter(
      comparison == "Tgfb1_Target_vs_NTC",
      group == "Egr2",
      abs(logFC) < 10
    )
  
  data_background <- 
    data_vis %>% 
    filter(!gene %in% genes)
  
  data_labels <- 
    data_vis %>% 
    filter(gene %in% genes)
  
  ggplot() +
    aes(logFC, -log10(p_adj)) +
    stat_binhex(
      data = data_background,
      aes(fill = "bg"),
      bins = 100,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c(bg = "lightgray")) +
    geom_point(
      data = data_labels,
      color = "blue",
      size = .05,
      shape = 20
    ) +
    geom_text_repel(
      data = data_labels,
      aes(label = gene),
      force_pull = 0.1,
      seed = 0,
      min.segment.length = 0,
      size = BASE_TEXT_SIZE_MM,
      segment.size = BASE_LINEWIDTH
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linewidth = BASE_LINEWIDTH,
      linetype = "dashed"
    ) +
    ylab("-log10 padj") +
    theme_pub() +
    theme(panel.grid = element_blank())
}

plot_volcano(labeled_genes)
ggsave_default("5f_volcano", type = "pdf", width = 50, height = 50)
# ggsave_default("5f_volcano_huge", type = "pdf", width = 400, height = 400)
