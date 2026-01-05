# -----------------------------------------------------------------------------
# Interactive Gene Expression Visualization Tool
# Collaborative sections clearly marked for easy editing by each team member
# -----------------------------------------------------------------------------

# Required libraries
library(shiny)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(plotly)
library(patchwork) # Used for combining plots
library(RColorBrewer)
library(data.table)
library(shinycssloaders)
library(htmlwidgets) # <-- ADD THIS LINE for saving HTML plots
library(colorspace)
library(Seurat) # for pseudobulking (https://satijalab.org/seurat/reference/aggregateexpression)
library(pheatmap)

# Helper functions for violin plots an box plots
get_y_coords <- function(i, n_cols, n_rows) {
  # Calculate which row this index is in (0-based)
  row <- (i - 1) %/% n_cols
  
  # Return y-coordinate based on row position
  # Top row (row 0) gets 1, bottom row gets 1/n_rows
  if (n_rows == 1) {
    return(rep(1, length(i)))
  }
  
  # Calculate y-coordinate: 1 - (row / (n_rows - 1)) * (1 - 1/n_rows)
  # This gives us the pattern: 1, 1, 0.67, 0.67, 0.33 for 5 items in 3 rows
  return(1 - (row / (n_rows - 1)) * (1 - 1/n_rows))
}

#' Get y-coordinate for subplot annotation (auto-calculate rows)
#' 
#' @param i Index of the subplot (1-based)
#' @param n_cols Number of columns in the subplot grid
#' @param total_items Total number of items (to calculate rows)
#' @return Y-coordinate for the annotation (0 to 1, where 1 is top)
get_y_coords_auto <- function(i, n_cols, total_items) {
  n_rows <- ceiling(total_items / n_cols)
  return(get_y_coords(i, n_cols, n_rows))
}

# Increase upload limit
options(shiny.maxRequestSize = 10000 * 1024^2)
Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)

# --- UI DEFINITION -----------------------------------------------------------
ui <- fluidPage(
  titlePanel("Interactive Gene Expression Visualization Tool - Sohum, Suki, Kim, and Zach"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_type", "What type of data will be used?",
                   choices = c("User-Uploaded" = "userload",
                               "Preloaded" ="preload"),
                   selected = "userload",
                   # inline = FALSE,
                   # width = NULL,
                   # choiceNames = NULL,
                   # choiceValues = NULL
      ),
      # fileInput("expr", "1. Expression Matrix (.tsv or .tsv.gz)", accept = c(".tsv", ".tsv.gz")),
      # fileInput("meta", "2. Metadata File (.tsv or .tsv.gz)",   accept = c(".tsv", ".tsv.gz")),
      
      # Show file inputs if user selects "User-Uploaded"
      conditionalPanel(
        condition = "input.data_type == 'userload'",
        fileInput("expr", "1. Expression Matrix (.tsv or .tsv.gz)", 
                  accept = c(".tsv", ".tsv.gz")),
        fileInput("meta", "2. Metadata File (.tsv or .tsv.gz)",   
                  accept = c(".tsv", ".tsv.gz"))
      ),
      
      # Show dropdown if user selects "Preloaded"
      conditionalPanel(
        condition = "input.data_type == 'preload'",
        selectInput("preload_choice", "Select preloaded dataset:",
                    choices = c("Option 1", "Option 2", "Option 3"))
      ),
      
      hr(),
      uiOutput("gene_ui"),      # Gene selector
      uiOutput("group_ui"),     # Group-by selector
      uiOutput("split_ui"),     # Split-by selector
      # checkboxInput("color_by_gene",
      #               "Feature Plot: Color by gene expression",
      #               value = FALSE),
      p(strong("Note:"), "For Violin/Box plots, multiple genes will be shown in a grid."),
      
      # --- NEW: Download Center UI ---
      hr(),
      # wellPanel(
      #   h4("Download Center"),
      #   p("Generate and download the combined UMAP feature plots below."),
      #   downloadButton("downloadPng", "Download as PNG", icon = icon("download")),
      #   downloadButton("downloadHtml", "Download as HTML", icon = icon("download"))
      # )
    ),
    mainPanel(
      tabsetPanel(
        # tabPanel("Feature Plot", withSpinner(plotlyOutput("featPlot", height = "900px"), type = 6)),
        # Separate the feature plot and split-plots
        tabPanel("Feature Plot", 
                 withSpinner(plotlyOutput("featPlot"), type = 6),
                 uiOutput("split_plot_ui")
        ),
        tabPanel("Bubble Plot",  withSpinner(uiOutput("dotPlot_ui"), type = 6)),
        tabPanel("Violin Plot",  withSpinner(uiOutput("vlnPlot_ui"), type = 6)),
        tabPanel("Heatmap",      plotlyOutput("heatmapPlot", height = "800px"), plotlyOutput("PBHeatmapPlot", height = "800px")),
        # tabPanel("Box Plot",     plotOutput("boxPlot",     height = "800px"))
        # Need a uiOutput for the boxplots that doesn't squeeze plots into a limited amount of space
        tabPanel("Box Plot",     uiOutput("boxPlotUI"))
      )
    )
  )
)

# --- SERVER LOGIC -----------------------------------------------------------
server <- function(input, output, session) {
  
  # 1) Load and merge expression + metadata -------------------------------
  plot_data <- reactive({
    # Instead of using as.matrix (which is slow), use data.table (which is way faster)!
    req(input$expr, input$meta)
    print("Checkpoint 1: plot_data")
    
    start_time <- Sys.time()
    
    # Transpose using data.table - this is the fastest method
    # keep.names creates a new column with the original column names
    # make.names uses the row names as new column names
    transposed_dt <- transpose(expr_dt_full(), keep.names = "cell", make.names = 'genes') 
    setkey(transposed_dt, cell)
    
    result <- transposed_dt[meta_dt_full(), nomatch = 0L]
    
    end_time <- Sys.time()
    print(paste("Time taken:", end_time - start_time))
    
    # Convert to data.frame if your downstream code requires it
    result <- as.data.frame(result)
  })
  
  pseudobulk_data = reactive({
    req(input$expr, input$meta, input$genes)
    shiny::validate(shiny::need(length(input$genes) > 1, "Must select more than one gene for pseudo-bulking"))
    
    expr = expr_dt_full()
    meta = meta_dt_full()
    
    # Clean column names
    setnames(expr, 1, "genes")
    setnames(meta, 1, "cell")
    expr[, genes := toupper(trimws(genes))]
    
    cat("Expression data:", nrow(expr), "genes x", ncol(expr)-1, "cells\n")
    cat("Metadata:", nrow(meta), "cells x", ncol(meta), "variables\n")
    
    # Convert to data.frame for easier handling
    # expr_df <- as.data.frame(expr)
    # meta_df <- as.data.frame(meta)
    
    # Why not filter the data down to just the genes you've selected? No need to process the entire matrix at once.
    expr_df <- expr[genes %in% input$genes]
    meta_df <- meta
    # browser()
    
    # Melt expression data to long format
    # Don't ude reshape2, use data.table
    expr_long <- data.table::melt(expr_df, id.vars = "genes", variable.name = "cell", value.name = "expression")
    expr_long$cell <- trimws(as.character(expr_long$cell))
    
    # Merge with metadata
    merged_data <- merge(expr_long, meta_df, by = "cell", all = FALSE)
    
    cat("Merged data:", nrow(merged_data), "observations\n")
    
    # Convert expression data to matrix format for Seurat
    expr_mat <- as.matrix(expr_df[, -1])  # Remove genes column
    rownames(expr_mat) <- expr_df$genes  # Set gene names as rownames
    colnames(expr_mat) <- colnames(expr_df)[-1]  # Preserve cell names
    
    cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "cells\n")
    
    # Choose grouping variable
    group_by_var <- input$group  # Change this as needed
    
    # Check if grouping variable exists and has valid values
    if (group_by_var %in% colnames(merged_data) && group_by_var != "None") {
      cat("Grouping by:", group_by_var, "\n")
      cat("Unique groups:", paste(unique(merged_data[[group_by_var]]), collapse = ", "), "\n")
      
      # Prepare metadata for Seurat - ensure cell names match
      meta_seurat <- meta_df
      rownames(meta_seurat) <- meta_seurat$cell
      
      # Verify cell names match between expression and metadata
      common_cells <- intersect(colnames(expr_mat), rownames(meta_seurat))
      cat("Common cells between expression and metadata:", length(common_cells), "\n")
      
      # Subset to common cells only
      expr_mat_subset <- expr_mat[, common_cells]
      meta_seurat_subset <- meta_seurat[common_cells, ]
      
      cat("Creating Seurat object...\n")
      # Create Seurat object
      seu <- CreateSeuratObject(
        counts = expr_mat_subset,
        meta.data = meta_seurat_subset
      )
      
      cat("Performing pseudobulking with AggregateExpression...\n")
      # Pseudobulking using Seurat's AggregateExpression
      pb <- AggregateExpression(
        seu,
        group.by = group_by_var,
        assays = "RNA",
        slot = "counts"
      )
      
      # Extract pseudobulk matrix and transpose
      pb_mat <- t(pb$RNA)
      
      # # Debug: Check what Seurat returned
      # cat("Raw Seurat output - pb$RNA dimensions:", nrow(pb$RNA), "x", ncol(pb$RNA), "\n")
      # cat("Raw Seurat output - pb$RNA rownames:", paste(head(rownames(pb$RNA)), collapse = ", "), "\n")
      # cat("Raw Seurat output - pb$RNA colnames:", paste(head(colnames(pb$RNA)), collapse = ", "), "\n")
      # 
      # # Ensure gene names are preserved
      # cat("After transpose - Gene names preserved:", length(colnames(pb_mat)), "genes\n")
      # cat("After transpose - Sample gene names:", paste(head(colnames(pb_mat)), collapse = ", "), "\n")
      # 
      # # Ensure group names are preserved
      # cat("After transpose - Group names preserved:", length(rownames(pb_mat)), "groups\n")
      # cat("After transpose - Group names:", paste(rownames(pb_mat), collapse = ", "), "\n")
      # 
      # cat("Pseudobulk matrix:", nrow(pb_mat), "groups x", ncol(pb_mat), "genes\n")
      
    } else {
      cat("No grouping variable available. Creating individual cell heatmap.\n")
      
      # No grouping: create matrix with individual cells
      # Select a subset of cells for visualization (to avoid too large matrices)
      max_cells <- 100  # Limit to 100 cells for performance
      if (ncol(expr_mat) > max_cells) {
        cat("Limiting to first", max_cells, "cells for visualization\n")
        expr_mat_subset <- expr_mat[, 1:max_cells]
      } else {
        expr_mat_subset <- expr_mat
      }
      
      # Transpose for individual cell visualization (cells as rows, genes as columns)
      pb_mat <- t(expr_mat_subset)
      
      # Ensure cell names are preserved
      cat("Cell names preserved:", length(rownames(pb_mat)), "cells\n")
      cat("Sample cell names:", paste(head(rownames(pb_mat)), collapse = ", "), "\n")
      
      # Ensure gene names are preserved
      cat("Gene names preserved:", length(colnames(pb_mat)), "genes\n")
      cat("Sample gene names:", paste(head(colnames(pb_mat)), collapse = ", "), "\n")
      
      # Update group_by_var for plotting
      group_by_var <- "Individual_Cells"
      
      cat("Individual cell matrix:", nrow(pb_mat), "cells x", ncol(pb_mat), "genes\n")
    }
    
    # print("pb_mat")
    # print(pb_mat)
    pb_mat
  })
  
  expr_dt_full <- reactive({
    req(input$expr)
    start_time = Sys.time()
    print("Checkpoint 2: expr_dt_full")
    dt <- fread(input$expr$datapath, encoding = "UTF-8", nThread = parallel::detectCores())
    setnames(dt, 1, "genes")
    dt[, genes := toupper(trimws(genes))]
    setkey(dt, genes)
    # print(colnames(dt))
    end_time = Sys.time()
    total_time = end_time - start_time
    print("Checkpoint 2 time")
    print(total_time)
    dt
  })
  
  meta_dt_full <- reactive({
    req(input$meta)
    print("Checkpoint 3: meta_dt_full")
    meta_dt <- fread(input$meta$datapath, data.table = TRUE)
    setnames(meta_dt, 1, "cell")
    meta_dt[, cell := trimws(as.character(cell))]
    setkey(meta_dt, cell)
    meta_dt
  })
  
  # 2) Dynamic UI for gene, group, split ----------------------------------
  output$gene_ui <- renderUI({
    req(expr_dt_full())
    print("Checkpoint 4: output$gene_ui")
    gene_choices <- unique(expr_dt_full()$gene)
    selectizeInput("genes", "Select Gene(s):",
                   choices  = NULL,
                   selected = NULL,
                   multiple = TRUE,
                   options = list(placeholder = "Type to search...", maxOptions = 100))
  })
  
  observeEvent(expr_dt_full(), {
    print("Checkpoint 5: observeEvent expr_dt_full")
    updateSelectizeInput(session, "genes",
                         choices = unique(expr_dt_full()$gene),
                         selected = NULL,
                         server = TRUE)
  })
  
  output$group_ui <- renderUI({
    req(meta_dt_full())
    print("Checkpoint 6: output$group_ui")
    meta_cols <- c("None", setdiff(names(meta_dt_full()), "cell"))
    
    # Remove UMAP_Xaxis and UMAP_Yaxis from meta_cols
    meta_cols <- setdiff(meta_cols, c("UMAP_Xaxis", "UMAP_Yaxis"))
    
    selectInput("group", "Group By (Categories/Colors):",
                choices  = meta_cols,
                selected = meta_cols[1])
  })
  output$split_ui <- renderUI({
    req(meta_dt_full())
    print("Checkpoint 7: output$split_ui")
    meta_cols <- setdiff(names(meta_dt_full()), c("cell", "UMAP_Xaxis", "UMAP_Yaxis"))
    selectInput("split", "Split By (Rows/Panels):",
                choices  = c("None", meta_cols),
                selected = "None")
  })
  
  # === FEATURE PLOT SECTION START (RESTRUCTURED FOR DOWNLOADS) ===
  
  # --- NEW: Reactive ggplot object for sharing ---
  feature_plotly_object <- reactive({
    
    req(plot_data(), "UMAP_Xaxis" %in% names(plot_data()), input$group)
    start_time = Sys.time()
    df <- plot_data()
    
    # --------------------------
    # Guard inputs
    # --------------------------
    req(input$group)
    req(input$split)
    split_col <- input$split
    group_col <- input$group
    
    # --------------------------
    # Prepare hover text
    # --------------------------
    # It should not matter if a group variable is selected or not
    if (!is.null(input$genes) && length(input$genes) > 0 && input$group != "None") {
      if (length(input$genes) == 1) {
        gene <- input$genes[1]
        summary_stats <- df %>%
          group_by(.data[[group_col]]) %>%
          summarise(
            avg_expr = mean(.data[[gene]], na.rm = TRUE),
            pct_expr = 100 * sum(.data[[gene]] > 0, na.rm = TRUE)/n(),
            .groups = "drop"
          )
        df <- left_join(df, summary_stats, by = group_col)
        df$hover_text <- paste0(
          "<b>Cell:</b> ", df$cell, "<br>",
          "<b>Group:</b> ", df[[group_col]], "<br>",
          "<b>Gene:</b> ", gene, "<br>",
          "<b>Expression (this cell):</b> ", round(df[[gene]],2), "<br>",
          "<b>Avg. Expr (group):</b> ", round(df$avg_expr,2), "<br>",
          "<b>% Expressing (group):</b> ", round(df$pct_expr,1), "%"
        )
      } else {
        df$ModuleScore <- rowMeans(df[, input$genes, drop = FALSE], na.rm = TRUE)
        summary_stats <- df %>%
          group_by(.data[[group_col]]) %>%
          summarise(avg_score = mean(ModuleScore, na.rm = TRUE), .groups = "drop")
        df <- left_join(df, summary_stats, by = group_col)
        df$hover_text <- paste0(
          "<b>Cell:</b> ", df$cell, "<br>",
          "<b>Group:</b> ", df[[group_col]], "<br>",
          "<b>Module Score (this cell):</b> ", round(df$ModuleScore,2), "<br>",
          "<b>Avg. Score (group):</b> ", round(df$avg_score,2)
        )
      }
    } else {
      # print("There was no group column")
      df$hover_text <- paste0("<b>Cell:</b> ", df$cell, "<br><b>Group:</b> ", df[[group_col]])
    }
    
    # --------------------------
    # Top panel: metadata UMAP
    # --------------------------
    # Color by selected group, if applicable
    if (input$group != "None"){
      meta_plotly <- plot_ly(
        data = df,
        x = ~UMAP_Xaxis,
        y = ~UMAP_Yaxis,
        color = ~.data[[group_col]], # Color by group
        text = ~hover_text,
        type = "scatter",
        mode = "markers",
        marker = list(size = 4, opacity = 0.8)
      ) %>%
        layout(
          title = "UMAP Colored by Metadata Group",
          legend = list(title = list(text = group_col))
        )
    } else {
      meta_plotly <- plot_ly(
        data = df,
        x = ~UMAP_Xaxis,
        y = ~UMAP_Yaxis,
        color = "red", # No groups to color by
        text = ~hover_text,
        type = "scatter",
        mode = "markers",
        marker = list(size = 4, opacity = 0.8)
      ) %>%
        layout(
          title = "UMAP Colored by Metadata Group",
          legend = list(title = list(text = group_col))
        )
    }
    
    # --------------------------
    # Bottom panel: split-feature UMAP
    # --------------------------
    split_plotly <- NULL
    if (!is.null(input$genes) && length(input$genes) > 0 && input$split != "None") {
      df_long <- df %>%
        select(UMAP_Xaxis, UMAP_Yaxis, cell, all_of(group_col), all_of(split_col), all_of(input$genes)) %>%
        pivot_longer(
          cols = all_of(input$genes),
          names_to = "gene",
          values_to = "expression"
        )
      
      # Background points
      df_bg <- df_long %>%
        filter(!is.na(.data[[split_col]]), !is.na(gene)) %>%
        distinct(UMAP_Xaxis, UMAP_Yaxis, .data[[split_col]], gene)
      
      # Hover text for split panel
      summary_split <- df_long %>%
        group_by(gene, .data[[split_col]]) %>%
        summarise(
          avg_expr_split = mean(expression, na.rm = TRUE),
          pct_expr_split = 100 * sum(expression > 0, na.rm = TRUE)/n(),
          .groups = "drop"
        )
      
      df_long <- left_join(df_long, summary_split, by = c("gene", split_col))
      df_long$hover_text_split <- paste0(
        "<b>Cell:</b> ", df_long$cell, "<br>",
        "<b>Group:</b> ", df_long[[group_col]], "<br>",
        "<b>", split_col, ":</b> ", df_long[[split_col]], "<br>",
        "<b>Gene:</b> ", df_long$gene, "<br>",
        "<b>Expression (this cell):</b> ", round(df_long$expression,2), "<br>",
        "<b>Avg. Expr (this panel):</b> ", round(df_long$avg_expr_split,2), "<br>",
        "<b>% Expressing (this panel):</b> ", round(df_long$pct_expr_split,1), "%"
      )
      
      facet_rows <- unique(df_long[[split_col]])
      facet_cols <- unique(df_long$gene)
      split_subplots <- list()
      
      # --------------------------
      # Create subplots
      # --------------------------
      for (g in facet_cols) {
        # If g was named "gene", bg_points and fg_points will not separate the dots by gene across plots. Don't make the iterating variable name the same as the column name that's being subsetted.
        print(g)
        row_subplots <- list()
        for (split_val in facet_rows) {
          print("split_val")
          # print(split_val)
          # bg_points = background “all cells” layer in light grey
          # fg_points = foreground “expressing cells” layer in color
          bg_points <- df_bg %>% filter(gene == g, .data[[split_col]] == split_val)
          fg_points <- df_long %>% filter(gene == g, .data[[split_col]] == split_val, expression > 0)
          
          print("fg_points")
          # print(fg_points)
          traces <- list()
          
          # This code creates an extra subplot that isn't needed
          # if (nrow(bg_points) > 0) {
          #   traces[[1]] <- plot_ly(
          #     data = bg_points,
          #     x = ~UMAP_Xaxis,
          #     y = ~UMAP_Yaxis,
          #     type = "scatter",
          #     mode = "markers",
          #     marker = list(size = 3, color = "grey90", opacity = 0.5),
          #     showlegend = FALSE,
          #     text = ~"",
          #     hoverinfo = "text"
          #   )
          # }
          
          # Interactive subplots.
          if (nrow(fg_points) > 0) {
            traces[[2]] <- plot_ly(
              data = fg_points,
              x = ~UMAP_Xaxis,
              y = ~UMAP_Yaxis,
              type = "scatter",
              mode = "markers",
              marker = list(
                size = 4,
                color = ~expression,
                colorscale = c("lightgrey","darkviolet"),
                showscale = (g == facet_cols[1] && split_val == facet_rows[1]), # only first plot
                colorbar = list(
                  thickness = 15,
                  len = 0.5,
                  y = 0.5,
                  title = list(text = "Expr")
                )
              ),
              text = ~hover_text_split,
              hoverinfo = "text",
              showlegend = FALSE,
              scaleratio = 1
            )
          }
          
          valid_traces <- traces[!sapply(traces, is.null)]
          if (length(valid_traces) > 0) {
            row_subplots[[split_val]] <- subplot(valid_traces, shareX = TRUE, shareY = TRUE)
          } else {
            row_subplots[[split_val]] <- NULL
          }
        }
        
        valid_rows <- row_subplots[!sapply(row_subplots, is.null)]
        if (length(valid_rows) > 0) {
          split_subplots[[g]] <- subplot(valid_rows, nrows = length(valid_rows), shareX = TRUE, shareY = TRUE)
        } else {
          split_subplots[[g]] <- NULL
        }
      }
      
      valid_cols <- split_subplots[!sapply(split_subplots, is.null)]
      if (length(valid_cols) > 0) {
        # split_plotly <- subplot(valid_cols, nrows = 1, shareX = FALSE, shareY = FALSE)
        
        # Keep the plots as square as possible (not squished)
        base_size <- 250
        split_plotly <- subplot(valid_cols, nrows = 1) %>%
          layout(
            autosize = FALSE,
            width  = base_size * length(facet_cols),
            height = base_size * length(facet_rows),
            xaxis = list(scaleanchor = "y"),
            yaxis = list(constrain = "domain", scaleratio = 1)
          )
        
        # --------------------------
        # Add facet labels
        # --------------------------
        annotations <- list()
        n_genes <- length(valid_cols)
        n_splits <- length(facet_rows)
        
        # Column labels (genes)
        for (i in seq_along(facet_cols)) {
          annotations <- c(annotations, list(
            list(x = (i - 0.5)/n_genes, y = 1.02, text = facet_cols[i],
                 xref = "paper", yref = "paper", showarrow = FALSE,
                 font = list(size = 12, face = "bold"))
          ))
        }
        
        # Row labels (split values)
        for (i in seq_along(facet_rows)) {
          annotations <- c(annotations, list(
            list(x = -0.02, y = 1 - (i - 0.5)/n_splits, text = facet_rows[i],
                 xref = "paper", yref = "paper", showarrow = FALSE,
                 font = list(size = 12, face = "bold"))
          ))
          # annotations <- c(annotations, list(
          #   list(
          #     x = -0.06,   # a bit outside the plots; adjust as needed
          #     y = 1 - (i - 0.5)/n_splits,
          #     text = facet_rows[i],
          #     xref = "paper", yref = "paper",
          #     showarrow = FALSE,
          #     textangle = -90,       # rotate text
          #     xanchor = "center",    # keep centered after rotation
          #     yanchor = "middle",
          #     font = list(size = 12, face = "bold")
          #   )
          # ))
        }
        
        split_plotly <- split_plotly %>% layout(annotations = annotations)
        
        # --------------------------
        # Combine with top panel
        # --------------------------
        combined_plot <- subplot(
          meta_plotly, split_plotly,
          nrows = 2, heights = c(0.3, 0.7),
          shareX = FALSE, shareY = FALSE,
          titleY = TRUE
        )
        
      } else {
        combined_plot <- meta_plotly
      }
      
    } else {
      combined_plot <- meta_plotly
    }
    
    end_time = Sys.time()
    time_taken = end_time - start_time
    print("Feature plot time")
    print(time_taken)
    
    # The feature plot and split-plots are separate elements now.
    return(list(meta_plotly, split_plotly))
  })
  
  # Render the plotly object for the UI
  output$featPlot <- renderPlotly({
    plot_gg <- feature_plotly_object()[[1]]
    print("Checkpoint 16: output$featPlot")
    req(plot_gg)
    ggplotly(plot_gg, tooltip = "text")
    # ggplotly(plot_gg)
    # return(plot_gg)
  })
  
  output$splitPlot <- renderPlotly({
    req(input$split)
    plot_gg <- feature_plotly_object()[[2]]
    print("Checkpoint 16A: output$splitPlot")
    req(plot_gg)
    ggplotly(plot_gg, tooltip = "text")
    # ggplotly(plot_gg)
    # return(plot_gg)
  })
  
  # Dynamically change the size of the panel based on how many split plots are made.
  # We don't want the plots to get "squished" horizontally.
  output$split_plot_ui <- renderUI({
    req(input$split) # this menu won't even appear until expression matrix and metadata are selected
    n_pixels_per_plot <- 200
    withSpinner(plotlyOutput("splitPlot", height = n_unique_splits() * n_pixels_per_plot), type = 6)    
  })
  
  
  # Calculate the number of rows in the bottom panel of the feature plot (split plots)
  n_unique_splits <- reactive({
    length(unique(plot_data()[[input$split]]))
  })
  
  # --- NEW: Download Handlers ---
  # output$downloadPng <- downloadHandler(
  #   filename = function() {
  #     paste0("feature_plot_", Sys.Date(), ".png")
  #   },
  #   content = function(file) {
  #     plot_to_save <- feature_plot_object()
  #     req(plot_to_save)
  #     # Use ggsave for high-quality export. Adjust dimensions as needed.
  #     ggsave(file, plot = plot_to_save, device = "png", width = 10, height = 12, dpi = 300, units = "in")
  #   }
  # )
  # 
  # output$downloadHtml <- downloadHandler(
  #   filename = function() {
  #     paste0("feature_plot_", Sys.Date(), ".html")
  #   },
  #   content = function(file) {
  #     plot_gg <- feature_plot_object()
  #     req(plot_gg)
  #     # Convert to plotly and then save as a self-contained HTML file
  #     plot_ly <- ggplotly(plot_gg, tooltip = "text")
  #     saveWidget(widget = plot_ly, file = file, selfcontained = TRUE)
  #   }
  # )
  
  # === FEATURE PLOT SECTION END ===
  
  # === BUBBLE PLOT SECTION START ===
  output$dotPlot_ui <- renderUI({
    
    if (is.null(input$split)) {
      n <- 1
    } else {
      n <- length(unique(dotPlot_data_prep()[[input$split]]))
    }
    
    cols <- 1
    rows <- ceiling(n/cols)
    
    height_px <- paste0(400*rows, "px")
    plotlyOutput("dotPlot", height = height_px, width = "100%")
  })
  
  dotPlot_data_prep = reactive({
    print("Approaching Checkpoint 17")
    print("Checkpoint 17: output$dotPlot")
    df <- plot_data()
    
    split_var <- if (input$split != "None") input$split else NULL
    group_var <- input$group
    
    dot_data <- df %>%
      select(all_of(c(input$genes, split_var, group_var))) %>%
      pivot_longer(cols = all_of(input$genes), names_to = "gene", values_to = "expression")
    
    if (!is.null(split_var)) {
      print("There is a split_var")
      summary_df <- dot_data %>%
        group_by(gene, .data[[split_var]], .add = TRUE) %>%
        { if (!is.null(group_var)) group_by(., .data[[group_var]], .add = TRUE) else . } %>%
        summarise(
          avg_expr = mean(expression, na.rm = TRUE),
          pct_expr = 100 * sum(expression > 0, na.rm = TRUE) / n(),
          .groups = "drop"
        )
    } else {
      # If there is no variable to split by, avoid an "object Object" error
      print("There is not a split_var")
      summary_df <- dot_data %>%
        group_by(gene, .data[[group_var]], .add = TRUE) %>%
        { if (!is.null(group_var)) group_by(., .data[[group_var]], .add = TRUE) else . } %>%
        summarise(
          avg_expr = mean(expression, na.rm = TRUE),
          pct_expr = 100 * sum(expression > 0, na.rm = TRUE) / n(),
          .groups = "drop"
        )
    }
    
    print("summary_df")
    print(summary_df)
    
    return(summary_df)
  })
  
  output$dotPlot <- renderPlotly({
    # if (input$split == "None" & !is.null(group_var)){
    #   split_var = group_var
    # }
    req(input$genes, input$group, plot_data())
    split_var <- if (input$split != "None") input$split else NULL
    group_var <- input$group
    
    start_time = Sys.time()
    # req(input$split != "None")
    summary_df = dotPlot_data_prep()
    # Hover text
    summary_df$hover_text <- paste0(
      "Gene: ", summary_df$gene, "<br>",
      if (!is.null(split_var)) paste0(split_var, ": ", summary_df[[split_var]], "<br>") else "",
      if (!is.null(group_var)) paste0(group_var, ": ", summary_df[[group_var]], "<br>") else "",
      "Avg. Expr: ", round(summary_df$avg_expr, 2), "<br>",
      "% Expr: ", round(summary_df$pct_expr, 1), "%"
    )
    
    print("summary_df")
    print(summary_df)
    
    ### bubble plot debug
    ### Initialize variables for creating subplots
    subplots <- list()
    genes <- unique(summary_df$gene)
    n_genes <- length(genes)
    annotations <- list()
    n_cols = 1
    
    # If there's 3 values or less in the group variable, bring the dots closer together.
    n_groups <- dplyr::n_distinct(summary_df[[group_var]])
    x_domain <- if (n_groups <= 3) c(0.20, 0.80) else c(0, 1)
    
    if (is.null(split_var)) {
      # No split_var included
      subplots[[1]] <- plot_ly(
        data = summary_df,
        x = ~.data[[group_var]], # group_var on x-axis
        y = ~gene, # genes on y-axis
        size = ~pct_expr,
        color = ~avg_expr,
        text = ~hover_text,
        hoverinfo = "text",
        type = "scatter",
        mode = "markers",
        marker = list(
          sizemode = "diameter",
          sizeref = 2,
          sizemin = 2,
          sizemax = 20,
          opacity = 0.9,
          line = list(width = 0.5, color = "white")
        ),
        colors = viridis::viridis_pal()(100)
      ) %>%
        layout(
          xaxis = list(
            title = input$group,
            type = "category",
            domain = x_domain
          ),
          yaxis = list(title = "Gene")
        )
    } else {
      # split_var included
      
      # Find the minimum and maximum avg_expr in the summary_df
      # to set the limits of the color bar legend.
      min_max <- range(summary_df$avg_expr)
      
      list_of_dfs <- as.data.table(summary_df) %>%
        split(by = split_var) # one graph for every split variable
      
      # This is for the plot labels
      n_split = length(list_of_dfs)
      
      # print("group_var")
      # print(group_var)
      # print("list_of_dfs")
      # print(list_of_dfs)
      # browser()
      
      # Color bar object. If we don't have this, the color bar will span the entire height of the plot.
      cb <- list(
        title = "avg_expr",
        lenmode = "fraction",
        len = 0.45,      # 45% of the subplot height
        y = 0.5,
        yanchor = "middle"
      )
      
      for (i in seq_along(list_of_dfs)) {
        subplots[[i]] <- plot_ly(
          data = list_of_dfs[[i]],
          # x = ~.data[[split_var]], # put the split_var on the x-axis
          # y = ~.data[[group_var]], # put the group_var on the y-axis
          x = ~.data[[group_var]], # put the group_var on the x-axis
          y = ~gene, # put the genes on the y-axis
          size = ~pct_expr,
          # color = ~avg_expr,
          text = ~hover_text,
          hoverinfo = "text",
          type = "scatter",
          mode = 'markers',
          marker = list(
            color = ~avg_expr,
            colorscale = 'Viridis',
            cmin = min_max[1],
            cmax = min_max[2],
            colorbar = cb,
            zauto = FALSE
          )
        ) %>%
          layout(
            xaxis = list(
              title = input$group,
              tickangle = 45,
              type = "category",
              domain = x_domain
            ),
            yaxis = list(title = "Gene")
          )
        
        annotations[[i]] <- list(
          x = 0.1,
          y = get_y_coords_auto(i, n_cols, n_split),
          text = unique(list_of_dfs[[i]][[split_var]]), # split value
          xref = "paper",
          yref = "paper", 
          xanchor = "center",
          yanchor = "bottom",
          showarrow = FALSE
        )
      }
      
    } #/ closes if is.null(split_var))
    
    end_time = Sys.time()
    time_taken = end_time - start_time
    print("Feature plot time")
    print(time_taken)
    
    # subplot(subplots, nrows = 1, shareY = TRUE, titleX = TRUE)
    if (length(subplots) > 1) {
      nrow <- ceiling(length(subplots) / n_cols)
      combined_plot <- subplot(subplots, nrows = nrow, shareX = TRUE, shareY = TRUE) %>%
        layout(annotations = annotations,
               title = "Gene Expression Bubble Plots",
               showlegend = TRUE,
               legend = list(orientation = "h", y = -0.1, x = 0.5, xanchor = "center")
        )
    } else {
      combined_plot <- subplots[[1]]
    }
  })
  
  # === BUBBLE PLOT SECTION END ===
  
  # === VIOLIN PLOT SECTION START (Enhanced) ===
  
  # Helper function for better colors
  get_colors <- function(n) {
    print("Checkpoint 18: get_colors")
    if (n <= 8) {
      brewer.pal(n, "Pastel1")
    } else if (n <= 12) {
      brewer.pal(n, "Set3")
    } else {
      qualitative_hcl(n, palette = "Pastel 1")
    }
  }
  
  output$vlnPlot_ui <- renderUI({
    req(input$genes)
    print("Checkpoint 19: output$vlnPlot_ui")
    n    <- length(input$genes)
    # cols <- if (n <= 2) n else 2
    cols <- 1
    rows <- ceiling(n/cols)
    height_px <- paste0(400*rows, "px")
    plotlyOutput("vlnPlot", height = height_px, width = "100%")
  })
  
  violinPlot_data = reactive({
    req(input$genes, input$group, expr_dt_full(), meta_dt_full())
    start_time = Sys.time()
    print("Checkpoint 20: output$vlnPlot")
    genes_selected <- toupper(trimws(input$genes))
    grp_col <- input$group
    spl_col <- input$split
    n_genes <- length(genes_selected)
    ncol_wrap <- if (n_genes <= 2) n_genes else 2
    
    withProgress(message = "Building violin plot(s)…", value = 0, {
      incProgress(0.25, detail = "Reading expression data...")
      dt_s <- expr_dt_full()[genes %in% genes_selected]
      
      if (nrow(dt_s) == 0) {
        showNotification("Selected genes not found in the expression matrix.", type = "error")
        return(NULL)
      }
      
      incProgress(0.1, detail = "Melting to long form")
      df_long <- melt(dt_s,
                      id.vars = "genes",
                      variable.name = "cell",
                      value.name = "expression",
                      variable.factor = FALSE)
      df_long[, cell := trimws(as.character(cell))]
      
      incProgress(0.1, detail = "Joining metadata")
      df <- merge(df_long, meta_dt_full(), by = "cell", all = FALSE)
      
      if (nrow(df) == 0) {
        showNotification("No matching cells between expression and metadata.", type = "error")
        return(NULL)
      }
      
      # Enhanced tooltip creation with metadata
      meta_cols <- setdiff(names(meta_dt_full()), "cell")
      df[, tooltip_text := {
        meta_info <- sapply(meta_cols, function(col) paste0(col, ": ", get(col)))
        paste(c(
          paste("Cell:", cell),
          paste0("Expression: ", round(expression, 3)),
          meta_info
        ), collapse = "<br>")
      }, by = seq_len(nrow(df))]
      
      incProgress(0.1, detail = "Filtering for densities")
      df[, group := df[[grp_col]]]
      if (spl_col != "None" && nzchar(spl_col)) {
        df[, split := df[[spl_col]]]
      }
      if (spl_col != "None" && nzchar(spl_col)) {
        cnts <- df[, .N, by = .(genes, group, split)]
        keep <- cnts[N >= 2, .(genes, group, split)]
        df_big <- merge(df, keep, by = c("genes", "group", "split"))
      } else {
        cnts <- df[, .N, by = .(genes, group)]
        keep <- cnts[N >= 2, .(genes, group)]
        df_big <- merge(df, keep, by = c("genes", "group"))
      }
      
      df_big[, group := factor(get(grp_col))]
      if (spl_col != "None" && nzchar(spl_col)) {
        df_big[, split := factor(get(spl_col))]
      } else {
        df_big[, split := NULL]
      }
      
      if (nrow(df_big) == 0) {
        showNotification("No groups have ≥2 cells after filtering. Try different settings.", type = "error")
        return(NULL)
      }
      end_time = Sys.time()
      time_taken = end_time - start_time
      print("Violin plot part 1 time")
      print(time_taken)
      # browser();
      return(df_big)
    })
  })
  
  output$vlnPlot <- renderPlotly({
    df_big = violinPlot_data()
    print("df_big")
    print(df_big)
    start_time = Sys.time()
    # Build plot
    genes <- unique(df_big$gene)
    print("genes")
    print(genes)
    n_genes <- length(genes)
    
    # Create subplots
    subplots <- list()
    annotations <- list()
    # n_cols <- 2
    n_cols <- 1
    for (i in seq_along(genes)) {
      print(paste0("genes[", i, "]"))
      print(genes[i])
      gene_data <- df_big[df_big$gene == genes[i], ]
      print("gene_data")
      print(gene_data)
      # browser()
      
      if (input$split != "None" && nzchar(input$split)) {
        # With split column
        gene_data$split = gene_data[[input$split]]
        p <- plot_ly(data = gene_data, showlegend = (i == 1)) %>%
          add_trace(
            x = ~group, y = ~expression, color = ~split,
            type = "violin",
            split = ~split,
            # width = 0.8,
            box = list(visible = TRUE),
            meanline = list(visible = TRUE),
            # fillcolor = ~split,
            line = list(color = "black"),
            opacity = 0.6,
            name = ~split,
            text = ~tooltip_text,
            hovertemplate = "%{text}<extra></extra>"
          ) %>%
          layout(
            violinmode = "group", # to split violins
            title = list(text = "", font = list(size = 12)),
            xaxis = list(title = "", tickangle = 45),
            yaxis = list(title = if(i == 1) "Normalized Expression" else ""),
            margin = list(l = 50, r = 20, t = 40, b = 50)
          )
        
        # %>%
        # add_trace(
        #   x = ~group, y = ~expression, color = ~split,
        #   type = "scatter", mode = "markers",
        #   marker = list(size = 3, opacity = 0.6),
        #   text = ~tooltip_text,
        #   hovertemplate = "%{text}<extra></extra>",
        #   showlegend = FALSE
        # )
      } else {
        # Without split column
        p <- plot_ly(data = gene_data, showlegend = (i == 1)) %>%
          add_trace(
            x = ~group, y = ~expression, color = ~group,
            type = "violin",
            box = list(visible = TRUE),
            meanline = list(visible = TRUE),
            # fillcolor = ~group,
            line = list(color = "black"),
            opacity = 0.6,
            name = ~group,
            text = ~tooltip_text,
            hovertemplate = "%{text}<extra></extra>"
          ) %>%
          layout(
            title = list(text = "", font = list(size = 12)),
            xaxis = list(title = "", tickangle = 45),
            yaxis = list(title = if(i == 1) "Normalized Expression" else ""),
            margin = list(l = 50, r = 20, t = 40, b = 50)
          )
      }
      
      subplots[[i]] <- p
      annotations[[i]] <- list(
        # x = ifelse(i %% 2 == 0, 0.8, 0.2),
        x = 0.5,
        y = get_y_coords_auto(i, n_cols, n_genes),
        text = genes[i],
        xref = "paper",
        yref = "paper",
        xanchor = "center",  
        yanchor = "bottom",  
        showarrow = FALSE
      )
    }
    
    # Create subplot with proper layout
    final_plot <- subplot(subplots, nrows = ceiling(n_genes / n_cols), 
                          shareY = FALSE, shareX = TRUE, 
                          titleX = FALSE, titleY = TRUE) %>%
      layout(
        annotations = annotations,
        margin = list(b = 150, t = 50),
        legend = list(orientation = "h", y = -0.15, x = 0.5, xanchor = "center"),
        showlegend = TRUE
      )
    
    end_time = Sys.time()
    time_taken = end_time - start_time
    print("Violin plot part 2 time")
    print(time_taken)
    return(final_plot)
  })
  # === VIOLIN PLOT SECTION END ===
  
  # === HEATMAP SECTION START ===
  # Function to min-max scale the values in the heatmap
  min_max_normalization <- function(x, min_value = -1, max_value = 1) {
    min_x <- min(x)
    max_x <- max(x)
    normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
    return(normalized_x)
  }
  
  output$heatmapPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    start_time = Sys.time()
    print("Checkpoint 21: output$heatmapPlot")
    
    # print(plot_data())
    
    df <- plot_data()
    # print(df)
    mat <- as.matrix(df[input$genes])
    # print(mat)
    print(class(mat))
    
    # Handle grouping
    if (input$group != "None") {
      grp <- df[[input$group]]
      ord <- order(grp)
      mat <- mat[ord, , drop = FALSE]
      grp_ordered <- grp[ord]
      
      avg_mat <- t(mat)
      
      group_levels <- unique(grp_ordered)
      group_bounds <- cumsum(table(grp_ordered))
      group_starts <- c(1, head(group_bounds + 1, -1))
      group_ends <- group_bounds
    } else {
      avg_mat <- t(mat)
      grp_ordered <- NULL
    }
    
    # Normalize (standardize genes)
    mat_z <- t(scale(t(avg_mat)))
    mat_z[is.na(mat_z)] <- 0
    
    # Min-max scale the heatmap rows
    normalized_data <- apply(mat_z, 1, min_max_normalization, min_value = -1, max_value = 1)
    print("Normalized data")
    norm_mat_z = t(normalized_data)
    # print(norm_mat_z)
    
    # print(max(mat_z))
    # print(min(mat_z))
    # print(pseudobulk_data())
    
    # Base heatmap
    p <- plot_ly(
      x = colnames(norm_mat_z),
      y = rownames(norm_mat_z),
      z = norm_mat_z,  # <- don't transpose here
      type = "heatmap",
      colors = colorRamp(c("violet", "black", "yellow")),
      showscale = TRUE
    )
    
    # If grouped, add labels + colored bars underneath
    if (input$group != "None") {
      palette <- RColorBrewer::brewer.pal(length(group_levels), "Set2")
      
      # Create colored bars for each group
      group_bars <- lapply(seq_along(group_levels), function(i) {
        list(
          type = "rect",
          x0 = group_starts[i],
          x1 = group_ends[i],
          y0 = -0.03, y1 = -0.015,
          xref = "x", yref = "paper",
          fillcolor = palette[i],
          line = list(color = palette[i])
        )
      })
      
      # Add text labels centered under each bar
      # Add text labels centered under each bar, rotated diagonally
      group_labels <- lapply(seq_along(group_levels), function(i) {
        list(
          x = mean(c(group_starts[i], group_ends[i])),
          y = -0.06,
          xref = "x", yref = "paper",
          text = group_levels[i],
          showarrow = FALSE,
          font = list(size = 12),
          textangle = 45,   # ← Rotate text diagonally
          xanchor = "center",
          yanchor = "top"
        )
      })
      
      p <- p %>% layout(
        title = list(
          text = "Expression Heatmap<br><sub>Cells grouped by variable</sub>",
          x = 0.5, xanchor = "center"
        ),
        xaxis = list(title = input$group, tickangle = 45, showticklabels = FALSE),
        yaxis = list(title = "Gene"),
        margin = list(l = 80, r = 20, b = 140, t = 80),
        shapes = group_bars,
        annotations = group_labels
      )
    } else {
      # No grouping: simple layout
      p <- p %>% layout(
        title = list(text = "Expression Heatmap", x = 0.5, xanchor = "center"),
        xaxis = list(title = "Cell", tickangle = 45),
        yaxis = list(title = "Gene"),
        margin = list(l = 80, r = 20, b = 100, t = 80)
      )
    }
    
    end_time = Sys.time()
    print("Heatmap time")
    print(end_time - start_time)
    
    p
  })
  
  # Pseudobulk heatmap
  output$PBHeatmapPlot <- renderPlotly({
    pb_mat = pseudobulk_data()
    selected_genes <- input$genes
    # available_genes <- intersect(selected_genes, colnames(pb_mat))
    group_by_var = input$group
    
    print("pb_mat")
    print(pb_mat)
    
    # Subset matrix for selected genes
    # heatmap_data <- pb_mat[, available_genes, drop = FALSE]
    
    # Actually, the data was subsetted in the pseudobulk_data() reactive
    heatmap_data <- pb_mat
    
    cat("Heatmap data:", nrow(heatmap_data), "groups x", ncol(heatmap_data), "genes\n")
    cat("Expression range:", round(min(heatmap_data), 3), "to", round(max(heatmap_data), 3), "\n")
    
    # Debug: Check matrix structure before melting
    cat("Matrix rownames (groups):", paste(rownames(heatmap_data), collapse = ", "), "\n")
    cat("Matrix colnames (genes):", paste(colnames(heatmap_data), collapse = ", "), "\n")
    
    # Ensure matrix has proper rownames and colnames
    if (is.null(rownames(heatmap_data)) || all(is.na(rownames(heatmap_data)))) {
      rownames(heatmap_data) <- paste0("Group_", 1:nrow(heatmap_data))
      cat("Added default rownames\n")
    }
    
    if (is.null(colnames(heatmap_data)) || all(is.na(colnames(heatmap_data)))) {
      colnames(heatmap_data) <- paste0("Gene_", 1:ncol(heatmap_data))
      cat("Added default colnames\n")
    }
    
    # Convert to long format for plotly
    # Create a data frame with proper factor levels to preserve names
    heatmap_df <- data.table(
      Group = rep(rownames(heatmap_data), each = ncol(heatmap_data)),
      Gene = rep(colnames(heatmap_data), times = nrow(heatmap_data)),
      Expression = as.vector(heatmap_data)
    )
    
    # Need to min-max scale the values in the pseudobulk heatmap (-1 to 1)
    # This code performs the normalization on a data frame instead of a matrix
    group_stats <- heatmap_df[, .(
      min_expr = min(Expression, na.rm = TRUE), 
      max_expr = max(Expression, na.rm = TRUE)
    ), by = Gene]
    
    min_value = -1; max_value = 1
    heatmap_with_stats_df = merge(heatmap_df, group_stats, by = "Gene")    
    heatmap_with_stats_df[min_expr != max_expr, normalized_expr := (Expression - min_expr) / (max_expr - min_expr) * (max_value - min_value) + min_value]
    
    # Ensure Group and Gene columns are character (not numeric)
    heatmap_with_stats_df$Group <- as.character(heatmap_with_stats_df$Group)
    heatmap_with_stats_df$Gene <- as.character(heatmap_with_stats_df$Gene)
    
    cat("Heatmap data frame dimensions:", nrow(heatmap_df), "observations\n")
    cat("Sample groups:", paste(unique(heatmap_df$Group), collapse = ", "), "\n")
    cat("Sample genes:", paste(unique(heatmap_df$Gene), collapse = ", "), "\n")
    
    # Create plotly heatmap
    if (group_by_var == "None") {
      # Individual cell heatmap
      plotly_heatmap <- plot_ly(
        data = heatmap_with_stats_df,
        x = ~Group,
        y = ~Gene,
        z = ~normalized_expr,
        type = "heatmap",
        colors = colorRamp(c("violet", "black", "yellow")),
        showscale = TRUE,
        colorbar = list(title = "Expression"),
        hovertemplate = paste(
          "<b>Cell:</b> %{y}<br>",
          "<b>Gene:</b> %{x}<br>",
          "<b>Expression:</b> %{z:.3f}<br>",
          "<extra></extra>"
        )
      ) %>%
        layout(
          title = list(
            text = "Individual Cell Expression Heatmap",
            x = 0.5,
            xanchor = "center"
          ),
          xaxis = list(
            # title = "Genes",
            tickangle = 45,
            showgrid = FALSE
          ),
          yaxis = list(
            title = "Genes",
            showgrid = FALSE
          ),
          margin = list(l = 80, r = 20, b = 100, t = 80)
        )
    } else {
      # Grouped pseudobulk heatmap
      plotly_heatmap <- plot_ly(
        data = heatmap_with_stats_df,
        x = ~Group,
        y = ~Gene,
        z = ~normalized_expr,
        type = "heatmap",
        colors = colorRamp(c("violet", "black", "yellow")),
        showscale = TRUE,
        colorbar = list(title = "Expression"),
        hovertemplate = paste(
          "<b>Group:</b> %{y}<br>",
          "<b>Gene:</b> %{x}<br>",
          "<b>Expression:</b> %{z:.3f}<br>",
          "<extra></extra>"
        )
      ) %>%
        layout(
          title = list(
            text = paste("Pseudobulk Heatmap"),
            x = 0.5,
            xanchor = "center"
          ),
          xaxis = list(
            title = "Groups",
            tickangle = 45,
            showgrid = FALSE
          ),
          yaxis = list(
            title = "Genes",
            showgrid = FALSE
          ),
          margin = list(l = 80, r = 20, b = 100, t = 80)
        )
    }
  })
  # === HEATMAP SECTION END ===
  
  # === BOX PLOT SECTION START ===
  # Need to change the height of the boxPlot plotOutput based on how many genes are selected
  output$boxPlotUI <- renderUI({
    if (is.null(input$genes) || length(input$genes) == 0) {
      return(NULL)   # no plot UI if nothing selected
    }
    nrows <- ceiling(length(input$genes) / 2)
    plotlyOutput("boxPlot", height = paste0(500 * nrows, "px"))
  })
  
  # This box plots takes pseudobulked data instead of the plot data
  # Actually, it can't. Pseudobulking computes single aggregated values for each group-variable and gene. I can't plot single data points in boxplots.
  output$boxPlot <- renderPlotly({
    req(input$genes, input$group, plot_data())
    req(input$group != "None")
    print("plot_data class")
    print(class(plot_data()))
    print(plot_data())
    
    # Actually, do not concert this to a data.frame. That's too time-consuming!
    # pseudobulk_df <- as.data.frame(pseudobulk_data())
    # df <- as.data.frame(plot_data())
    
    # print("pseudobulk as data frame")
    # print(class(pseudobulk_df))
    # print(pseudobulk_df)
    
    start_time = Sys.time()
    # For this boxplot, we need to pseudobulk the data first.
    # For pseudobulking there is a seurat function that does this well or you can make your own. Essentially it takes the individual cell ID columns form the expression matrix and the cell IDs as rows from the metadata file and averages the counts per biological replicate sample
    
    # df <- df %>%
    #   pivot_longer(cols = all_of(input$genes), names_to="gene", values_to="expr")
    
    # This code should be written in data.table, not dplyr.
    if (!is.null(input$split) && input$split != "None") {
      gene_df = plot_data() %>% select(all_of(input$genes), input$split, input$group) %>% as.data.table()  
    } else {
      gene_df = plot_data() %>% select(all_of(input$genes), input$group) %>% as.data.table()  
    }
    
    df <- data.table::melt(
      gene_df, # Only melt the genes that are selected. That will speed up runtime.
      measure.vars = input$genes,
      variable.name = "gene",
      value.name = "expr"
    )
    
    print("new df")
    
    print(df)
    
    df$hover_text <- paste0(
      "Gene: ", df$gene, "<br>",
      "Group: ", df[[input$group]], "<br>",
      if (input$split != "None") paste0("Split: ", df[[input$split]], "<br>") else "",
      "Expression: ", round(df$expr, 3)
    )
    
    # Create subplots for each gene
    subplots <- list()
    genes <- unique(df$gene)
    
    n_genes <- length(genes)
    n_cols <- 2
    annotations <- list()
    for (i in seq_along(genes)) {
      gene_data <- df[df$gene == genes[i], ]
      
      if (input$split != "None") {
        # With split variable
        p <- plot_ly(
          data = gene_data,
          x = ~.data[[input$group]],
          y = ~expr,
          color = ~.data[[input$split]],
          type = "box",
          text = ~hover_text,
          hoverinfo = "text",
          showlegend = (i == 1)
        ) %>%
          layout(
            boxmode = 'group', # this splits the boxes based on split-by feature
            xaxis = list(title = input$group, tickangle = 45),
            yaxis = list(title = genes[i])
          )
      } else {
        # Without split variable
        p <- plot_ly(
          data = gene_data,
          x = ~.data[[input$group]],
          y = ~expr,
          color = ~.data[[input$group]],
          type = "box",
          text = ~hover_text,
          hoverinfo = "text",
          showlegend = FALSE
        ) %>%
          layout(
            boxmode = 'group', # this splits the boxes based on split-by feature
            xaxis = list(title = input$group, tickangle = 45),
            yaxis = list(title = genes[i])
          )
      }
      
      subplots[[i]] <- p
      annotations[[i]] <- list(
        x = ifelse(i %% 2 == 0, 0.8, 0.2),
        y = get_y_coords_auto(i, n_cols, n_genes),
        text = genes[i],
        xref = "paper",
        yref = "paper",
        xanchor = "center",  
        yanchor = "bottom",  
        showarrow = FALSE
      )
    }
    
    
    # Combine subplots
    if (length(subplots) > 1) {
      ncol <- 2
      nrow <- ceiling(length(subplots) / ncol)
      combined_plot <- subplot(subplots, nrows = nrow, shareX = TRUE, shareY = FALSE) %>%
        layout(annotations = annotations,
               title = "Gene Expression Box Plots",
               showlegend = TRUE,
               legend = list(orientation = "h", y = -0.1, x = 0.5, xanchor = "center")
        )
    } else {
      combined_plot <- subplots[[1]]
    }
    
    end_time = Sys.time()
    time_taken = end_time - start_time
    print("Boxplot time")
    print(time_taken)
    return(combined_plot)
  })
  # === BOX PLOT SECTION END ===
}

# Launch the app
shinyApp(ui, server)