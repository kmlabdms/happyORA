library(tidyverse)
library(ComplexHeatmap)


# Process signatures ------------------------------------------------------
path_sheet = "./data-raw/Transcription Signature Gene Lists (2).xlsx"
sheets = readxl::excel_sheets(path_sheet)
gene_list = map(sheets, ~readxl::read_xlsx(path_sheet, sheet = .x, col_names = F) %>%
                    rename(gene_id = `...1`)) %>% set_names(sheets)


# RUNNING -----------------------------------------------------------------
# example of query_input
comb1 <- bind_rows(gene_list[1:2], .id = "cluster")
comb2 <- bind_rows(gene_list[3:6], .id = "cluster")
comb3 <- bind_rows(gene_list[7:8], .id = "cluster")
col1 = c("#D95F02","#1B9E77")
col2 = c("#7296B9", "#98276C", "#e29357", "#d76a9a")
col3 = c("#F8766D","#00BFC4")


# happyORA output a plot and a table for overlapped genes.
g1 = happyORA::happyORA(query_input = comb1, colors = col1)
g2 = happyORA(query_input = comb2, colors = col2)
g3 = happyORA(comb3, colors = col3)
# g2 plot is too tall so we split it
g2_nohallmark = happyORA::happyORA(comb2, col2,
                                   group_labels = c(unique(comb2$cluster), "XXX"),
                                   remove_genesets_group = "Hallmark",
                                   visual_query_as_row = T)
g2_nohallmark =  happyORA(comb2, col2, keep_genesets_group  = "Hallmark")

# Output overlapped genes
list(comb1 = g1$tbl, comb2 = g2$tbl, comb3 = g3$tbl) %>%
  writexl::write_xlsx("output/overlapped_genes_3_comb.xlsx")


cols_us= c("#7296B9", "#98276C", "#e29357", "#d76a9a","#F8766D","#00BFC4")
gene_group = read_tsv("/Users/qixu/Github/mycmiecases/gene_to_gsea.tsv")
g1 = happyORA(query_input = gene_group %>% filter(str_detect(cluster, "WGS")),
              remove_genesets_group = c("Hallmark"),
              colors = unlist(map(cols_us,~rep(.x, 1)))
              )
count(gene_group, cluster)


cbio = cBioPortalData::cBioPortal()
gene_fmi =
  cBioPortalData::getGenePanel(cbio, "FoundationOne") %>%
  pull(hugoGeneSymbol)

g2 = happyORA(query_input = gene_group %>% filter(str_detect(cluster, "Panel")),
              remove_genesets_group = c("Amplicon"),
              background = gene_fmi,
              colors = unlist(map(cols_us,~rep(.x, 1)))
)
count(gene_group, cluster)



## CCLE cell line driver genes

cl_driver_path = "data-raw//gdsc/cell_line_genomic_feature_clean.rds"
cl_driver_raw = read_rds(cl_driver_path) %>%
  select(cl_id, is_mutated, gain_loss, genes)


cl_mapping = read_csv("data-raw//gdsc/sample_info.csv") %>%
  select(CCLE_Name, cl_id = COSMICID) %>%
  filter(!grepl("Merged", CCLE_Name)) %>%
  dplyr::rename(sample_id = CCLE_Name)

cl_driver <-
  cl_driver_raw %>%
  left_join(cl_mapping) %>%
  drop_na(sample_id) %>%
  separate_rows(genes) %>%
  select(sample_id, genes)

q_samples = c("TGBC11TKB_STOMACH", "SNU216_STOMACH", "2313287_STOMACH", "SNU1_STOMACH", "SNU520_STOMACH", "FU97_STOMACH", "GSS_STOMACH", "MKN7_STOMACH", "MKN74_STOMACH", "NUGC2_STOMACH", "ECC12_STOMACH", "MKN45_STOMACH", "HGC27_STOMACH", "IM95_STOMACH", "KATOIII_STOMACH")



q_signature =
  cl_driver %>%
  filter(sample_id %in% q_samples) %>%
  rename(cluster = sample_id, gene_id = genes)


ora_res = happyORA::happyORA(query_input = q_signature,
                             group_labels = q_samples,
                             remove_genesets_group = "Hallmark",
                             visual_query_as_row = T)
(gsea_plot = ora_res$plot)







