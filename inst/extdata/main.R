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
g2 = happyORA(comb2, colors = col2)
g3 = happyORA(comb3, colors = col3)
# g2 plot is too tall so we split it
g2_nohallmark = happyORA(comb2, col2, remove_genesets_group = "Hallmark")
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








