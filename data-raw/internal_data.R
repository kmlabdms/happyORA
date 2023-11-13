gene_gw <- readRDS("data-raw/gsea/gene_gw_all.rds")

geneset_type_all <-
  readxl::read_xlsx("data-raw/gsea/ref_genesets_230418.xlsx",
                    sheet = "ref_genesets_230327_all")

usethis::use_data(geneset_type_all, gene_gw, internal = T, overwrite = T)
