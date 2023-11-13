#' SimpleGSEA
#'
#' @param query_input two column tibble with "gene_id" and "cluster" as group
#' @param genesets_list  lists of gene sets
#' @param background can be genome-wide. if not provide,
#' bg will be total genes of genesets_list
#' @param p_cut 0.05 by default
#' @param cols  colors in order
#'
#' @return
#' @export
#'
#' @examples
SimpleGSEA <- function(query_input, genesets_list, background = NULL, visual = T,
                       p_cut = 0.05, cols = NULL) {

  #geneset_type <- read_csv(paste0(dir, "data/gsea_data/geneset_type.csv"))

  gene_signautes <- genesets_list
  # if (is.null(background)) background <- readRDS("./data/gsea/gene_gw_all.rds")

  # Main Func
  gses_result <- runCancerGSEA(
    query_input,
    gene_signautes,
    gene_panel = background
  )

  # Visual
  if (visual) {
    df <- gses_result %>% gather(group, value, -geneset)
    pbase <-
      ggplot(df, aes(x = group, y = str_wrap(geneset, 50))) +
      geom_point(size = 6, color = "lightgrey") +
      geom_point(data = df[df$value < 0.05, ], aes(color = group), size = 6) +
      labs(  x = "", y = "" )
    pformat <- pbase +
      theme(axis.text = element_text(size = 10, face = "bold")) +
      theme(legend.position = "none")
    if (!is.null(cols)){
      pformat = pformat + scale_color_manual(values = cols)
    }
  } else {
    pformat = NULL
  }



  # formating
  gses_genes <-
    gses_result %>%
    gather("type", "pvalue", -geneset) %>%
    arrange(pvalue) %>%
    filter(pvalue <= p_cut) %>%
    mutate(overlapped = map2_chr(
      .$geneset, .$type, .f = OverlapGene,
      # extract parameter
      all_signature = gene_signautes,
      all_query = query_input
    )) %>%
    arrange(geneset, type) %>%
    mutate(pvalue = round(pvalue, 3)) %>%
    mutate(geneset = fct_infreq(geneset) ) %>%
    arrange(geneset, type)

  output <- list(gses_genes = gses_genes, gses_plot = pformat)
  return(output)

}




