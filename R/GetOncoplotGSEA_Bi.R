#' GetOncoplotGSEA_Bi
#'
#' @param tbl_overlap
#' @param tbl_spl
#' @param group_colors
#' @param group_labels
#' @param name
#' @param visual_query_as_row
#' @param ...
#'
#' @return
#' @importFrom ComplexHeatmap alter_graphic oncoPrint
#' @export
#'
#' @examples
GetOncoplotGSEA_Bi <- function(tbl_overlap,
                            tbl_spl = NULL,
                            group_colors = c("tomato1", "limegreen", "blue"),
                            group_labels = NULL,
                            name = "p<0.05",
                            visual_query_as_row = F,
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
    full[missing_group] <- ""
  }

  mat <- column_to_rownames(full, "geneset")
  if (!visual_query_as_row) {

    mat <- mat[, group_labels]

    row_split <- tibble(geneset = rownames(mat)) %>% left_join(tbl_spl)
    ht <- oncoPrint(mat,
                    alter_fun = alt_fun, show_pct = F,
                    column_order = group_labels,
                    width = unit(5, "mm") * ncol(mat),
                    height = unit(5, "mm") * nrow(mat),
                    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                    column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                    #row_gap = unit(5, "mm"),
                    border = T, column_names_rot = 90,
                    top_annotation =NULL, right_annotation = NULL,
                    row_split = row_split[["type"]],  row_title_rot = 0,
                    heatmap_legend_param = list(title = name),
                    alter_fun_is_vectorized = FALSE,
                    ...)
  } else {
    mat <- t(mat) # sample as rows
    mat <- mat[group_labels, ]
    col_split <- tibble(geneset = colnames(mat)) %>% left_join(tbl_spl) %>%
      mutate(type = str_wrap(type, 5))
    ht <- oncoPrint(mat,
                    alter_fun = alt_fun, show_pct = F,
                    width = unit(6, "mm") * ncol(mat),
                    height = unit(5, "mm") * nrow(mat),
                    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                    column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),

                    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                    border = T,
                    row_order  = group_labels,
                    top_annotation =NULL, right_annotation = NULL,
                    column_split = col_split[["type"]],  column_title_rot = 0,
                    show_column_names = T,
                    heatmap_legend_param = list(title = name),
                    alter_fun_is_vectorized = FALSE)
  }


  return(ht)
}
