library(ComplexHeatmap)

GetOncoplotGSEA <- function(tbl_overlap,
                            tbl_spl = NULL,
                            group_colors = c("tomato1", "limegreen", "blue"),
                            group_labels = NULL,
                            name = "p<0.05",
                            ...){


  alt_fun = map(group_colors, ~alter_graphic("rect", fill = .x))
  names(alt_fun) = group_labels
  alt_fun = c(alt_fun, background = alter_graphic("rect", fill = "grey90"))


  tbl_overlap$geneset <- str_remove_all(tbl_overlap$geneset, "HALLMARK_")
  tbl_spl$geneset <- str_remove_all(tbl_spl$geneset, "HALLMARK_")



  full <- tbl_overlap %>%
    select(geneset, type) %>%
    mutate(value = as.character(type)) %>%
    spread(type, value, fill = "")

  missing_group <- setdiff(group_labels, colnames(full))
  if (length(missing_group) > 0 ) {
    full[[missing_group]] <- ""
  }

  mat <- column_to_rownames(full, "geneset")

  # mat <- tbl_overlap %>% select(geneset, type) %>%
  #   mutate(value = type) %>%
  #   spread(type, value, fill = "") %>%
  #   column_to_rownames("geneset")

  if (!is.null(tbl_spl)) {
    row_split <- tibble(geneset = rownames(mat)) %>% left_join(tbl_spl)
    ht <- oncoPrint(mat,
                    alter_fun = alt_fun, show_pct = F,
                    column_order = group_labels,
                    width = unit(5, "mm") * ncol(mat),
                    height = unit(5, "mm") * nrow(mat),
                    #row_gap = unit(5, "mm"),
                    border = T, column_names_rot = 90,
                    top_annotation =NULL, right_annotation = NULL,
                    row_split = row_split[["type"]],  row_title_rot = 0,
                    heatmap_legend_param = list(title = name),
                    alter_fun_is_vectorized = FALSE,
                    ...)
  } else {
    ht <- oncoPrint(mat, top_annotation =NULL,
                    right_annotation = NULL,
                    row_gap = unit(5, "mm"),
                    width = unit(5, "mm") * ncol(mat),
                    height = unit(5, "mm") * nrow(mat),
                    show_pct = F,
                    ...)
  }

  #p <- cheatmap2ggplot(ht)

  return(ht)
}





if (F) {
  source("./function/gsea/SimpleGSEA.R")
  # Genesets
  geneset_type_all <- readxl::read_xlsx("data/gsea/ref_genesets_230418.xlsx",sheet = "ref_genesets_230327_all")
  genesets <-
    geneset_type_all %>%
    separate_rows(genes) %>%
    split(.$geneset) %>%
    map(~pull(.x, genes))
  #genesets <- readRDS(paste0("./data/gsea/ref_genesets_230413.rds"))
  gw <- readRDS("./data/gsea/gene_gw_all.rds")
  bk <- read_csv("data/FM_panel/FM_gene_panel_tidy.csv")[["gene"]]

  selected <- c("Cancer_related",
                "Gene Targets",
                "Pathway Targets",
                "Other Targets")
  group_list <- geneset_type_all %>% split(.,.$type) %>% map(~pull(.x, geneset))
  geneset_type <- filter(geneset_type_all, type %in% selected)
  selector <- unlist(group_list[selected])
  genesets_selected <- genesets[selector]

  # Query

  query_input <- read_csv("data/shiny_input_colorectal_example.csv") %>%
      select(cluster = type, gene_id = gene)
  query_input <- bind_rows(query_input, mutate(query_input,cluster = "Total"))

  # tts_sig <- read_csv("/Users/qixu/Github/GenA/shiny-app/app/data/TTS-clustering/tts_sig_level_1_10C.csv") %>%
  #   select(cluster = sample_id, gene_id = gene)
  # query_input <- tts_sig

  # RUNNING
  gsea_result <- SimpleGSEA(query_input, genesets_selected,
                            background = gw, p_cut = 0.05)
  tbl_overlap <- gsea_result$gses_genes
  tbl_overlap$type <- factor(tbl_overlap$type, levels = unique(query_input$cluster))
  GetOncoplotGSEA(tbl_overlap, geneset_type)

  GetOncoplotGSEA_general(tbl_overlap, geneset_type)
  clipr::write_clip(tbl_overlap %>% arrange(type))

}


