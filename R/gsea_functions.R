library(patchwork)
dir <- "./"

# base GSEA function ----
do_hyper <- function(ids, universe, gene_sets) {
  retval <- as.data.frame(t(sapply(gene_sets, function(gene_set) {
    x <- sum(ids %in% gene_set)
    n <- length(ids)
    X <- sum(universe %in% gene_set)
    N <- length(universe)
    if (x > 0 && X > 0) {
      pval <- phyper(x - 1, X, N - X, n, lower.tail = FALSE)
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = (x / n) / (X / N),
        "pval" = pval
      )
    } else {
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = 0,
        "pval" = 1
      )
    }
    retvec
  })))
  retval$Name <- names(gene_sets)
  retval <- retval[order(retval$pval),]
  retval$qval <- p.adjust(retval$pval, method = "fdr")
  retval
}


# Build GSEA function ----
do_hyper_ma <-
  function(sig, sig_universe, gene_sets) {
    hyper_output <- do_hyper(
      ids = sig,
      universe = sig_universe,
      gene_sets = gene_sets
    )
    output <- hyper_output %>%
      select(pval) %>%
      as_tibble(rownames = "geneset")
    output
  }


# Run GSEA function and clean the data ----


runCancerGSEA <-
  function(signature, gene_sets, gene_panel = NULL) {
    
    # 3 Input query genes ----
    query <- arrange(signature, cluster)
    gene_list_four <- lapply(
      split(query, query$cluster),
      function(table) {
        pull(table, gene_id)
      }
    )

    # 4 Perform GSEA  ----
    output_ma <-
      lapply(
        X = gene_list_four, FUN = do_hyper_ma,
        sig_universe = gene_panel,
        gene_sets = gene_sets
      ) %>%
      plyr::join_all(by = "geneset", type = "inner") %>%
      purrr::set_names("geneset", names(gene_list_four)) %>%
      as_tibble()
    output_ma
  }


# Run ggplot for GSEA visualization ----
gsea_ggplot <- function(
  gses_result, pvalue_cut, project_name = "GSEA", geneset_type, 
  level_order = NULL, color_order = NULL, show_p = F, maintitle = NULL
  ) {
  
  
  
  gsea_pvalue <- gses_result %>% inner_join(geneset_type)
  data <- gsea_pvalue %>%
    gather("group", "pvalue", -geneset, -type) %>%
    mutate(pvalue = ifelse(pvalue <= pvalue_cut, pvalue, NA))
  
  if (!is.null(level_order)) {
    data$group <- factor(data$group , levels = level_order)
    color_names <- level_order
  }
  
  if ((!is.null(level_order)) & any(str_detect(colnames(gsea_pvalue), "Total"))) {
    total_title = str_extract(colnames(gsea_pvalue), "VUS") %>% na.omit() %>% unique()
    data$group <- factor(data$group , 
                         levels = c(level_order, total_title))
    color_names <- c(level_order, total_title)
  }
  
  if (is.null(level_order)) {color_names <- unique(data$group)}
  
  if (is.null(color_order)) {
    additional_colors <- 
      try(c( "red", "purple", "orange","blue",  "limegreen" )[1:length(color_names)])
  } else {
    additional_colors <- color_order
  }
  
  names(additional_colors) <- color_names
  
  data_rename <- data %>%
    mutate(geneset = geneset %>%
             str_replace_all("_", " ") %>%
             str_to_title()) %>%
    arrange(group, pvalue) %>%
    mutate(geneset = replace(geneset, geneset == "Tsg", "TSG")) %>% 
    mutate(geneset = factor(geneset, levels = rev(unique(geneset))))
  yaxis_label <- "Gene sets"
  
  ## 5.1 cancer biology  -----
  if (FALSE) {
    data_cancer <-
      filter(data_rename, type %in% c("Hallmark")) %>%
      filter(pvalue <= pvalue_cut) %>%
      bind_rows(filter(data_rename, type %in% c("Cancer_related")))
  } 
 
  geneset_to_show <- 
    data_rename %>% 
    filter(pvalue < pvalue_cut | geneset %in% c("TSG", "Oncogene")) 
    
  
  data_cancer <- 
    data_rename %>% 
    filter(geneset %in% unique(geneset_to_show$geneset)) %>% 
    filter(type %in% c("Hallmark", "Cancer_related")) %>%
    mutate(pvalue <- ifelse(pvalue < pvalue_cut, pvalue, NA))
  

  
  gsea_dot1 <-
    ggplot(
      data_cancer,
      aes(
        x = group,
        y = str_wrap(geneset, 50),
        label = round(pvalue, 2)
      )
    ) +
    geom_point(size = 6, color = "lightgrey") +
    geom_point(data = na.omit(data_cancer), aes(color = group), size = 6) +
    scale_color_manual(values = c("Total" = "blue", "Main" = "red", "VUS" = "limegreen"))+
    #scale_color_manual(values = c("blue", "red", "limegreen", "orange", "purple"))+
    
    theme_set(theme_bw()) +
    theme(
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.text.x = element_text(
        size = 10,
        angle = 45, vjust = 1, hjust = 1
      ),
      axis.text.y = element_text(face = "bold", size = 8),
      title = element_text(face = "bold", size = 10)
    ) +
    theme(axis.ticks = element_blank()) +
    theme( axis.text.x = element_blank()) +
     theme(legend.position = "none")+
    
    facet_wrap("type", scales = "free_y", nrow = 1) +
    labs(
      title = glue::glue("Cancer Biology"),
      #subtitle = project_name,
      x = "", y = ""
    )
  #gsea_dot1
  if (!any(str_detect(data_cancer$group, "VUS"))) {
    gsea_dot1 <- gsea_dot1 + 
      scale_color_manual(values = additional_colors)
  }
  if (show_p) {
    gsea_dot1 <- gsea_dot1 +
      geom_text(color = "white", fontface = "bold", size = 2) 
  }
  

  
  ## 5.2 treatment   -----
  data_cancer <-
    filter(data_rename, !type %in% c("Hallmark", "Cancer_related", "Amplicon")) %>%
    mutate(geneset = factor(geneset, levels = rev(unique(geneset)))) %>%
    identity()
  
  gsea_dot2 <-
    ggplot(
      data_cancer,
      aes(x = group, y = str_wrap(geneset, 35), label = round(pvalue, 2))) +
    geom_point(size = 6, color = "lightgrey") +
    geom_point(data = na.omit(data_cancer), aes(color = group), size = 6) +
    scale_color_manual(values = c("Total" = "blue", "Main" = "red", "VUS" = "limegreen")) +
    
     theme_set(theme_bw()) +
    theme(
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.text.x = element_text(
        size = 8,
        angle = 45, vjust = 1, hjust = 1
      ),
      axis.text.y = element_text(face = "bold", size = 8),
      title = element_text(face = "bold", size = 10)
    ) +
    theme( axis.text.x = element_blank()) +
    theme(axis.ticks = element_blank()) +
    facet_wrap("type", scales = "free", nrow = 1) +
    #theme(legend.position = "none")+
    labs(
      title = glue::glue("Treatment related"),
      x = "", y = ""
    )
  #gsea_dot2
  if (!any(str_detect(data_cancer$group, "VUS"))) {
    gsea_dot2 <- gsea_dot2 + 
      scale_color_manual(values = additional_colors, name = maintitle)
  }
  if (show_p) {
    gsea_dot2 <- gsea_dot2 +
      geom_text(color = "white", fontface = "bold", size = 2)
  }

  
  data_amplicon <-
    filter(data_rename, type %in% c("Amplicon")) %>%
    identity() %>%
    mutate(group = factor(group, levels = rev(unique(group))))
  data_amplicon$geneset <-
    factor(data_amplicon$geneset,
           levels = unique(data_amplicon$geneset) %>% str_sort(numeric = T)
    )
  
  
  gsea_dot3 <-
    ggplot(
      data_amplicon,
      aes(
        x = group,
        y = geneset,
        label = round(pvalue, 2)
      )
    ) +
    geom_point(size = 6, color = "lightgrey") +
    geom_point(data = na.omit(data_amplicon), aes(color = group), size = 6) +
    scale_color_manual(values = c("Total" = "blue", "Main" = "red", "VUS" = "limegreen"))+
    theme_set(theme_bw()) +
    theme(
      axis.text = element_text(size = 10, face = "bold", color = "black"),
      axis.text.x = element_text(
        size = 8, 
        angle = 45, vjust = 1, hjust = 1
      ),
      axis.text.y = element_text(face = "bold", size = 8),
      
      title = element_text(face = "bold", size = 10)
    ) +
    theme(legend.position = "none")+
    
    theme(axis.ticks = element_blank()) +
    coord_flip() +
    facet_wrap("type", scales = "free", nrow = 1) +
    labs(y = "",  x = "")
  #gsea_dot3
  #gsea_dot2
  if (!any(str_detect(data_cancer$group, "VUS"))) {
    gsea_dot3 <- gsea_dot3 + 
      scale_color_manual(values = additional_colors)
  }
  
  if (show_p) {
    gsea_dot3 <- gsea_dot3 +
      geom_text(color = "white", fontface = "bold", size = 2)
  }
  
  
  
  p2 <- gsea_dot1 | (gsea_dot2 / gsea_dot3) +
    plot_layout(heights = c(3.5, 1), widths = c(1, 3))
  
  layout <- "
              AABBB
              AABBB
              AABBB
              AACCC
              "
  p3 <-
    gsea_dot1 + gsea_dot2 + gsea_dot3 +
    plot_layout(design = layout) +
    plot_annotation(title = maintitle) 
  
  return(p3)
}



# extrav specific genes ----
OverlapGene <- function(signature, query, all_signature, all_query) {
  gene_query <- filter(all_query, cluster %in% query) %>% pull(gene_id)
  gene_sig <- all_signature[[signature]]
  output <- paste(intersect(gene_query, gene_sig), sep = "/",collapse = "/")
  return(output)
}


#main function
mainGSEA <- function(query_genes, background_genes, pvaluecut = 0.2,  
                     add_total = TRUE,
                     level_order = NULL, color_order = NULL, show_p = FALSE,
                     maintitle = NULL
                     ) {
  
  geneset_type <- read_csv(paste0(dir, "data/gsea_data/geneset_type.csv"))
  gene_signautes <- readRDS(paste0(dir, "data/gsea_data/ma_list_211116.rds"))
  
  
  if (add_total == TRUE) {
    query_input <- 
      query_genes %>% 
      mutate(cluster = "Total") %>%  
      bind_rows(query_genes)
  } else {
    query_input <- query_genes
  }
  
  gses_result <- runCancerGSEA(query_input, gene_signautes,gene_panel = background_genes)
  gses_plot <- gsea_ggplot(gses_result, pvalue_cut=pvaluecut, 
                           geneset_type = geneset_type,
                           level_order = level_order, color_order = color_order,
                           show_p = show_p, maintitle = maintitle
                           )
  
  gses_genes <- 
    gses_result %>%  
    gather("type", "pvalue", -geneset) %>% 
    arrange(pvalue) %>% filter(pvalue < pvaluecut) %>% 
    mutate(overlapped = map2_chr(
      .$geneset, .$type, .f = OverlapGene,
      # extract parameter
      all_signature = gene_signautes, 
      all_query = query_input
    )) %>%
    arrange(geneset, type) %>% 
    mutate(pvalue = round(pvalue, 3))
  
  output <- list(gses_genes = gses_genes, gses_plot = gses_plot)
  return(output)
  
}



# Testing
test = FALSE
if (test == TRUE) {
  background_genes <- read_csv("data/FM_panel/FM_gene_panel_tidy.csv")[["gene"]]
  query_genes <-
    read_csv("data/shiny_input_colorectal_example.csv") %>%
    select(gene_id = gene, cluster = type)
  test_output <- mainGSEA(
    query_genes = query_genes,
    pvalue = 0.2, 
    background_genes = background_genes,
    add_total = F
    )
  test_output$gses_plot
}


if (test == TRUE) {
  background_genes <- readRDS("~/Dropbox/Kowalski/src/gene_gw_all.rds")
  query_genes <-
    read_csv("/Users/qixu/Dropbox/Kowalski/CRC_Subtyping/data/gene_signature_four.csv") %>%
    select(gene_id = gene, cluster = signature)
  test_output <- mainGSEA(
    query_genes,
    pvalue = 0.05, 
    background_genes,
    add_total = F
  )
  test_output$gses_plot
}

if (test == TRUE) {
  background_genes <- readRDS("~/Dropbox/Kowalski/src/gene_gw_all.rds")
  query_genes <- read_csv("/Users/qixu/Library/CloudStorage/Dropbox/Kowalski/CRC_Subtyping/data/table_of_four_sigature.csv") %>% 
    filter(signature == "Primary") %>% 
    select(gene_id = gene, cluster = cluster) 
  test_output <- mainGSEA(
    query_genes,
    pvalue = 0.05, 
    background_genes,
    add_total = F,
    maintitle = "PRimart"
  )
  test_output$gses_plot
  test_output$gses_genes %>% view()
}







