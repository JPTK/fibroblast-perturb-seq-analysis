library(SingleCellExperiment)
library(tidyverse)
library(ggrepel)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

col_data <-
  read_rds("data_generated/rna.rds") %>% 
  colData() %>% 
  as_tibble()



# Analyze -----------------------------------------------------------------

plot_ntc_depletion <- function() {
  plot_data <- 
    col_data %>% 
    summarise(
      .by = cell_type,
      n = n(),
      n_ntc = sum(target == "non-targeting"),
      n_ko = n - n_ntc
    )
  
  ggplot(plot_data, aes(n_ko + 1, n_ntc + 1)) +
    geom_point(
      aes(color = cell_type),
      size = 0.5
    ) +
    geom_abline(linewidth = BASE_LINEWIDTH) +
    scale_x_log10(
      "Number of cells with KO guide + 1",
      limits = c(1, NA)
    ) +
    scale_y_log10(
      "Number of cells with NTC guide + 1",
      limits = c(1, NA)
    ) +
    scale_color_manual(
      name = "Cell type",
      values = CELL_TYPE_COLORS
    ) +
    coord_fixed() +
    geom_text_repel(
      data = plot_data %>% filter(n_ko / 10 > n_ntc),
      aes(label = cell_type),
      size = BASE_TEXT_SIZE_MM
    ) +
    theme_pub() +
    theme(
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(2, "mm"),
      panel.grid = element_blank()
    )
}

plot_ntc_depletion()
ggsave_default("Sx_ntc_depletion", type = "pdf", width = 70, height = 50)
