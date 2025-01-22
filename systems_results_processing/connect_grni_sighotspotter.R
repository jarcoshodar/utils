library(dplyr)
library(igraph)

# A small helper to remove trailing "...2", "...21", etc. 
# from duplicated gene names that cause mismatches
clean_label <- function(x) {
  sub("\\.\\.\\.[0-9]+$", "", x)
}

#drop row names from a data frame also to avoid auto-labeled duplicates
drop_rownames <- function(df) {
  rownames(df) <- NULL
  df
}

parse_rds_to_list <- function(rds_path) {
  rds_data <- readRDS(rds_path)
  cell_type_names <- names(rds_data)
  
  gather_nodes_edges <- function(vis_list) {
    nd <- data.frame()
    ed <- data.frame()
    if (!is.null(vis_list) && length(vis_list) > 0) {
      for (hotspot_name in names(vis_list)) {
        hotspot <- vis_list[[hotspot_name]]
        if (!is.null(hotspot$nodes)) {
          tmp_nodes <- drop_rownames(hotspot$nodes)
          nd <- bind_rows(nd, tmp_nodes)
        }
        if (!is.null(hotspot$edges)) {
          tmp_edges <- drop_rownames(hotspot$edges)
          ed <- bind_rows(ed, tmp_edges)
        }
      }
    }
    list(nodes = nd, edges = ed)
  }
  
  parsed_list <- lapply(cell_type_names, function(ct) {
    ct_data <- rds_data[[ct]]
    
    # Combine from vis_net_A and vis_net_I
    a_res <- gather_nodes_edges(ct_data$vis_net_A)
    i_res <- gather_nodes_edges(ct_data$vis_net_I)
    
    all_nodes <- bind_rows(a_res$nodes, i_res$nodes) %>% distinct()
    all_edges <- bind_rows(a_res$edges, i_res$edges) %>% distinct()
    
    # Clean up trailing "...2", etc.
    if ("id" %in% names(all_nodes)) {
      all_nodes$id <- clean_label(all_nodes$id)
    }
    if ("from" %in% names(all_edges)) {
      all_edges$from <- clean_label(all_edges$from)
      all_edges$to   <- clean_label(all_edges$to)
    }
    
    #if the sighotspoter edges have a column named "Effect" or "effect" that is numeric,
    #rename it to "rds_effect" to avoid conflict with gurobi's textual "effect".
    if ("Effect" %in% names(all_edges)) {
      names(all_edges)[names(all_edges) == "Effect"] <- "rds_effect"
    }
    if ("effect" %in% names(all_edges)) {
      # If this is numeric or -1/1, rename it
      names(all_edges)[names(all_edges) == "effect"] <- "rds_effect"
    }
    
    # Also carry over the cell_type
    all_nodes$cell_type <- ct
    all_edges$cell_type <- ct
    
    list(nodes = all_nodes, edges = all_edges)
  })
  
  names(parsed_list) <- cell_type_names
  parsed_list
}


parse_csv_data <- function(csv_path) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  required_cols <- c("TF", "effect", "target", "Score")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("CSV is missing required columns: ",
         paste(missing_cols, collapse=", "))
  }
  
  #Each row => an edge from=TF, to=target
  #rename "effect" -> "csv_effect" to avoid conflict with sighotspotter
  edges <- df %>%
    rename(csv_effect = effect) %>%
    mutate(
      from = TF,
      to   = target,
      source = "CSV"
    ) %>%
    select(from, to, csv_effect, Score, source)
  
  #the node set is the union of all from/to
  nodes <- data.frame(
    id = unique(c(edges$from, edges$to)),
    source = "CSV",
    stringsAsFactors = FALSE
  )
  
  list(nodes = nodes, edges = edges)
}

build_combined_network <- function(rds_nodes, rds_edges, csv_nodes, csv_edges) {
  #Standardize effect notation
  standardize_effect <- function(df) {
    #for sighotspotter
    if ("rds_effect" %in% names(df)) {
      df <- df %>% 
        mutate(effect = case_when(
          rds_effect == 1 ~ "Activation",
          rds_effect == -1 ~ "Inhibition",
          TRUE ~ NA_character_
        ))
    }
    
    #for grnopt
    if ("csv_effect" %in% names(df)) {
      df <- df %>% 
        mutate(effect = case_when(
          is.na(effect) & csv_effect == "Activation" ~ "Activation",
          is.na(effect) & csv_effect == "Inhibition" ~ "Inhibition",
          TRUE ~ effect
        ))
    }
    
    df
  }
  
  #sighotspotter edges
  if (!("source" %in% names(rds_edges))) {
    rds_edges$source <- "RDS"
  }
  if (!("Score" %in% names(rds_edges))) {
    rds_edges$Score <- "no_data"  #"no_data" instead of NA
  }
  
  #grnopt edges
  if (!("cell_type" %in% names(csv_edges))) {
    csv_edges$cell_type <- NA
  }
  
  #combine and make similar
  all_edges <- bind_rows(
    rds_edges %>% 
      select(from, to, rds_effect, source, cell_type) %>%
      mutate(csv_effect = NA_character_),
    csv_edges %>% 
      select(from, to, csv_effect, Score, source) %>%
      mutate(rds_effect = NA_real_,
             Score = as.character(Score))  #Score to character for consistency
  ) %>%
    standardize_effect() %>%
    select(from, to, effect, Score, source, cell_type) %>%
    distinct()
  
  ##TODO: ensure metadata is properly kept, for now
  all_nodes <- data.frame(
    id = unique(c(rds_nodes$id, csv_nodes$id)),
    stringsAsFactors = FALSE
  )
  #in igraph, convert after return if you need dataframe
  g <- graph_from_data_frame(
    d = all_edges,
    directed = TRUE,
    vertices = all_nodes
  )
  
  #for plotting
  E(g)$color <- ifelse(E(g)$effect == "Activation", "green", "red")
  
  g
}
read_and_build_network <- function(rds_path, csv_path, cell_type) {
  
  #sighotspotter should have named elements based on cell type
  rds_list <- parse_rds_to_list(rds_path)
  
  #parse grni 
  csv_list <- parse_csv_data(csv_path)
  
  #act on mistmatch between network and sighotspotter
  if (! cell_type %in% names(rds_list)) {
    stop("Cell type not found in RDS: ", cell_type)
  }
  
  #get sighotspotter data for this celltype
  rds_nodes <- rds_list[[cell_type]]$nodes
  rds_edges <- rds_list[[cell_type]]$edges
  
  g <- build_combined_network(rds_nodes, rds_edges,
                              csv_list$nodes, csv_list$edges)
  g
}

test_graph <- read_and_build_network(
  rds_path = "~/Documents/Final_Mohammed/sigFINAL/Exercise2_results3_Cerebellum_sighotspotter_hotspots_filtered_elements_all.rds",
  csv_path = "~/Documents/Final_Mohammed/genenetworks4/Exercise2_results_Integrated_Cerebellum_GRNI5NeurpOE_vs__NeurpOC.csv",
  cell_type = "NeurpOE_vs_NeurpOC"
)

# Inspect the graph:
test_graph
as_data_frame(test_graph, "vertices")  #see the node set, probably only this interacts with intercom but check all
usable_result <- as_data_frame(test_graph, "edges")     # see columns like from, to, rds_effect, csv_effect, Score, ...
