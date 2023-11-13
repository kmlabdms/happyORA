#' happyORA
#'
#' @param query_input a two-column table with column name as "cluster" and "gene_id"
#' @param colors
#' @param remove_genesets_group a vector or a string. group names of gene sets not to show
#' "Amplicon", "Cancer_related", "CHIP", "DNA_Damage_Repair", "Epigenetic",
#' "Gene Targets", "Hallmark", "Other Targets", "Pathway Targets"
#' @param keep_genesets_group select a group of gene sets only
#' @param background if null, genome-wide genes are used
#' @param p_cut default 0.05
#'
#' @return
#' @export
#'
#' @examples
happyORA <- function(query_input,
                         colors = c("red", "blue", "limegreen"),
                         remove_genesets_group = NULL,
                         keep_genesets_group = NULL,
                         background = NULL,
                         p_cut = 0.05,
                     group_labels = NULL,
                     visual_query_as_row = F
                         ){

  # for testing only
  # remove_genesets_group = NULL
  # keep_genesets_group = NULL
  # background = NULL
  # p_cut = 0.05

  if (is.null(background)) {
    background <- gene_gw
  }

  if (length(unique(group_labels)) > length(colors)) {
    colors = rep(colors[1], length(group_labels))
  }

  # database Prep -----------------------------------------------------------
  # geneset_type_all <-
  #   readxl::read_xlsx("data/gsea/ref_genesets_230418.xlsx",sheet = "ref_genesets_230327_all")
  group_list <- geneset_type_all %>% split(.,.$type) %>% map(~pull(.x, geneset))
  genesets <-
    geneset_type_all %>%
    separate_rows(genes) %>%
    split(.$geneset) %>%
    map(~pull(.x, genes))


  genesets_group_selected = names(group_list)
  if (!is.null(remove_genesets_group)) {
    genesets_group_selected <-
      genesets_group_selected[!genesets_group_selected %in% remove_genesets_group]
  }
  if (!is.null(keep_genesets_group)) {
    genesets_group_selected = keep_genesets_group
  }

  selector <- unlist(group_list[genesets_group_selected])
  genesets_selected <- genesets[selector]
  geneset_type <- filter(geneset_type_all, type %in% genesets_group_selected)


  # ORA
  gsea_result <- SimpleGSEA(query_input, genesets_selected,
                            background = background, p_cut = p_cut)
  tbl_overlap <- gsea_result$gses_genes
  tbl_overlap$type <- factor(tbl_overlap$type,
                             levels = str_sort(unique(query_input$cluster), numeric = T))
  if (is.null(group_labels)) {
    group_labels = str_sort(unique(query_input$cluster), numeric = T)
  }
  ht = GetOncoplotGSEA_Bi(
    tbl_overlap = tbl_overlap,
    tbl_spl =  geneset_type,
    group_colors = colors,
    group_labels =  group_labels,
    show_column_names= T,
    visual_query_as_row = visual_query_as_row,
    name = sprintf("p<%s", p_cut)
                       )

  return(list(plot = ht, tbl = tbl_overlap))
}











