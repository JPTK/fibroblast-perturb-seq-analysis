library(fgsea)
library(tidyverse)
library(fs)
source("scripts_r/utils.R")



# Load data ---------------------------------------------------------------

dge <- 
  dir_ls("data_generated/", regex = "DEG_") %>% 
  read_csv(id = "file") %>% 
  separate_wider_regex(file, c(".+DEG_", comparison = ".*", "\\.csv")) %>% 
  select(comparison, group, gene = names,
         logFC = logfoldchanges, p_adj = pvals_adj)

gene_map <-
  read_tsv("metadata/biomart_human_mouse_20210727.tsv") %>% 
  select(human_gene = `Human gene name`, mouse_gene = `Gene name`) %>% 
  distinct()



# Perform enrichment ------------------------------------------------------

enrichr_databases <- c(
  # pathways
  "KEGG_2019_Mouse",
  "MSigDB_Hallmark_2020",
  "NCI-Nature_2016",
  "Reactome_2022",
  "WikiPathways_2019_Mouse",
  
  # transcription factor targets
  "ChEA_2022",
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "ENCODE_TF_ChIP-seq_2015",
  "TRANSFAC_and_JASPAR_PWMs",
  "TRRUST_Transcription_Factors_2019"
)

# download Enrichr databases in a format that can be used by fgsea:
# assemble a named list (names = databases)
#   of named lists (names = enrichment terms)
#   of character vectors (all genes associated with the respective term).
enrichr_genesets <-
  enrichr_databases %>% 
  set_names() %>% 
  map(\(db) {
    info("Downloading {db}")
    url <- paste0(
      "https://maayanlab.cloud/Enrichr/geneSetLibrary",
      "?mode=text&libraryName=",
      db
    )
    db <- read_lines(url)
    m <- str_match(db, "(.+?)\\t\\t(.+)")
    terms <- m[, 2]
    genes <- m[, 3] %>% str_split("\\t")
    genes %>% 
      map(stringi::stri_remove_empty) %>% 
      set_names(terms)
  })

# remove human gene sets from TRRUST database
enrichr_genesets$TRRUST_Transcription_Factors_2019 <-
  enrichr_genesets$TRRUST_Transcription_Factors_2019 %>% 
  magrittr::extract(imap_lgl(., ~str_detect(.y, "mouse"))) %>% 
  set_names(str_extract, "\\w+")

# convert human to mouse genes
enrichr_genesets <- 
  enrichr_genesets %>% 
  modify(\(db) {
    modify(db, \(gene_set) {
      gene_map %>% 
        filter(human_gene %in% gene_set) %>% 
        pull(mouse_gene) %>% 
        unique()   
    })
  })



# perform gene set enrichment analysis;
# returns dataframe, comprising columns "db", "comparison", "group",
# as well as all columns in the result of fgsea().
run_gsea <- function(comparison, group, db) {
  ranked_genes <-
    dge %>%
    filter(comparison == {{comparison}}, group == {{group}}) %>%
    arrange(desc(logFC)) %>% 
    select(gene, logFC) %>%
    deframe()
  ranked_genes <- ranked_genes[!is.na(ranked_genes)]
  
  # info("GSEA of comparison {comparison}, group {group}, ",
  #      "db {db} ({length(ranked_genes)} genes)")
  
  fgsea(
    enrichr_genesets[[db]],
    ranked_genes,
    eps = 0
  ) %>%
    as_tibble() %>%
    mutate(
      db = {{db}},
      comparison = {{comparison}},
      group = {{group}},
      .before = 1
    )
}

# run_gsea("Il1b_Target_vs_NTC", "Arid2", "KEGG_2019_Mouse")  # for testing


gsea_results <-
  expand_grid(
    comparison = unique(dge$comparison),
    group = unique(dge$group),
    db = names(enrichr_genesets)
  ) %>% 
  # head(5) %>%
  pmap(run_gsea, .progress = TRUE) %>% 
  list_rbind()



# Save results ------------------------------------------------------------

write_rds(dge, "data_generated/dge.rds")
write_rds(gsea_results, "data_generated/gsea.rds")
