library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyFiles)
library(shinycssloaders)
library(tidyverse)
library(plotly)
library(rlang)
library(viridis)
library(DT)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(writexl)


library(scico)

ui <- navbarPage("LupaCell",
                 theme = NULL, 
                 header = tags$head(
                   tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
                 ),
                 # Home
                 tabPanel("Home",
                          fluidPage(
                            titlePanel("Welcome to LupaCell"),
                            p("PPC"),
                            p("Use the tabs above to explore datasets, heatmaps, spatial plots, and more."),
                            hr(),
                            fluidRow(
                              column(6,
                                     actionButton("save_config", "Save Current Data Settings", class = "btn-success")
                              ),
                              column(6,
                                     fileInput("load_config", "Load Saved Settings", accept = ".rds")
                              )
                            )
                          )
                 ),
                 
                 # Upload Data
                 # tabPanel("Upload Data",
                 #          fluidPage(
                 #            fileInput("user_file", "Upload CSV, TSV, or RDS"),
                 #            verbatimTextOutput("file_info")
                 #          )
                 # ),
                 
                 tabPanel("Upload Data",
                          fluidPage(
                            fluidRow(
                              column(6,
                                     shinyFilesButton("user_file", "Choose CSV, TSV, or RDS File", "Browse", multiple = FALSE)
                              ),
                              column(6,
                                     shinyFilesButton("meta_file", "Load Metadata File", "Browse Metadata", multiple = FALSE)
                              )
                            ),
                            
                            h4("File Info"),
                            withSpinner(verbatimTextOutput("file_info"), type = 4)
                          )
                 ),
                 
                 
                 ## gene tab
                 tabPanel("Gene",
                          fluidPage(
                            h3("Gene Annotation Tool"),
                            actionButton("annotate_genes", "Annotate Gene Symbols", class = "btn-primary"),
                            br(),
                            actionButton("convert_to_symbol", "Convert to Gene Symbols", class = "btn-warning"),
                            br(), br(),
                            withSpinner(DTOutput("gene_annotation_table"), type = 4),
                            hr(),
                            h3("Highly Variable Gene Selector"),
                            numericInput("top_n_hvg", "Number of top variable genes:", value = 2000, min = 10, step = 100),
                            fluidRow(
                              column(6, actionButton("run_hvg", "Select HVGs", class = "btn-primary")),
                              column(6, checkboxInput("overwrite_data", "Overwrite with HVGs", value = FALSE)),
                              column(6, checkboxInput("overwrite_log_norm", "Overwrite with LogNorm", value = FALSE))
                            ),
                            downloadButton("download_hvg_list", "Download HVG List"),
                            br(), br(),
                            withSpinner(DTOutput("hvg_table"), type = 4),
                            br(),
                            withSpinner(plotlyOutput("hvg_plot"), type = 4)
                          )
                 ),
                 
                 
                 
                 
                 # Umap viewer
                 tabPanel("Plot Viewer",
                          fluidPage(
                            h3("Interactive Viewer"),
                            fluidRow(
                              column(4,
                                     uiOutput("umap_x_ui"),
                                     uiOutput("umap_y_ui"),
                                     uiOutput("umap_color_ui"),
                                     radioButtons("var_type", "Variable Type:",
                                                  choices = c("Categorical", "Numerical"),
                                                  inline = TRUE),
                                     actionButton("plot_umap", "Generate Plot", class = "btn-primary"),
                                     actionButton("save_umap_to_dashboard", "Add  Plot to Dashboard", class = "btn-success")
                              ),
                              column(8,
                                     plotlyOutput("umap_plot")
                              )
                            )
                          )
                 ),
                 
                 
                 tabPanel("Heatmap",
                          fluidPage(
                            h3("Heatmap Generator"),
                            pickerInput("genes_heatmap", "Select Genes for Heatmap:", choices = NULL, 
                                        multiple = TRUE, options = list(`live-search` = TRUE)),
                            selectInput("group_by_heatmap", "Group by Metadata Variable:", choices = NULL),
                            actionButton("generate_heatmap", "Generate Heatmap", class = "btn-primary"),
                            downloadButton("download_heatmap_matrix", "Download Scaled Matrix"),
                            br(), br(),
                            withSpinner(plotOutput("heatmap_plot", height = "800px"), type = 4)
                          )
                 ),
                 
                 tabPanel("Expression timeline",
                          fluidPage(
                            h3("Pseudotime Heatmap"),
                            selectInput("pseudo_var", "Select Pseudotime Variable:", choices = NULL),
                            sliderInput("timeline_bins", "Number of Pseudotime Bins:", min = 5, max = 300, value = 50, step = 5),
                            pickerInput("timeline_genes", "Select Genes:", choices = NULL, 
                                        multiple = TRUE, options = list(`live-search` = TRUE)),
                            actionButton("plot_timeline", "Plot Pseudotime Heatmap", class = "btn-primary"),
                            br(), br(),
                            withSpinner(plotOutput("timeline_heatmap", height = "800px"), type = 4)
                          )
                 ),
                 
                 tabPanel("Tissue Plot",
                          fluidPage(
                            h3("Spatial Gene Expression"),
                            selectInput("tissue_x", "X Axis (Coordinate):", choices = NULL),
                            selectInput("tissue_y", "Y Axis (Coordinate):", choices = NULL),
                            selectInput("tissue_gene", "Select Gene:", choices = NULL),
                            actionButton("save_tissue_to_dashboard", "Add Plot to Dashboard", class = "btn-success"),  # <- new button
                            withSpinner(plotlyOutput("tissue_plot"), type = 4)
                          )
                 ),
                 
                 tabPanel("Dashboard",
                          fluidPage(
                            h3("Saved Timeline Plot"),
                            plotOutput("saved_timeline_plot"),
                            h3("HVG Plot"),
                            plotlyOutput("saved_hvg_plot"),  # <- Add this line
                            uiOutput("dashboard_umaps"),
                            uiOutput("dashboard_tissue")
                          )
                 ),
                 
                 
                 
                 
                 
                 
                 absolutePanel(
                   top = 10, right = 20, width = "auto", draggable = FALSE,
                   style = "z-index:9999; background-color: rgba(255,255,255,0.8); padding: 5px 10px; border-radius: 5px; font-size: 14px; color: #444;",
                   textOutput("file_name_label")
                 )
                 
)












server <- function(input, output, session) {
  
  roots <- c(home = normalizePath("~"))
  
  uploaded_data <- reactiveVal(NULL)   # this way you can easily update the data by data(newdata). works really well
  metadata <- reactiveVal(NULL) 
  
  
  file_path = reactiveVal(NULL)
  
  file_name = reactiveVal(NULL)
  
  
  meta_path = reactiveVal(NULL)
  meta_file_name= reactiveVal(NULL)
  
  
  # for gene hgihlight
  selected_hvg_gene <- reactiveVal(NULL)
  
  
  
  ## saving plots 
  
  hvg_plot = reactiveVal(NULL)
  
  gene_expression  = reactiveVal(NULL)
  
  timeline_gene_expression = reactiveVal(NULL)
  

  umap_dashboard_plots <- reactiveVal(list())
  
  tissue_dashboard_plots <- reactiveVal(list())
  
  last_tissue_plot <- NULL
  ## track files
  config <- reactiveValues(
    count_path = NULL,
    meta_path = NULL,
    hvg_genes = NULL
  )
  
  # for above 
  observeEvent(selected_path(), {
    config$count_path <- selected_path()
  })
  
  observeEvent(metadata_path(), {
    config$meta_path <- metadata_path()
  })
  
  observeEvent(hvg_result(), {
    config$hvg_genes <- hvg_result()$top_genes
  })
  
  
  # save button
  observeEvent(input$save_config, {
    showModal(modalDialog("Saving configuration...", easyClose = TRUE))
    saveRDS(list(
      count_path = file_path(),
      file_name= file_name(),
      meta_path = meta_path(),
      hvg_genes = config$hvg_genes
    ), file = "lupacell_config.rds")
    showNotification("Settings saved to 'lupacell_config.rds'", type = "message")
  })
  
  
  
  observeEvent(input$load_config, {
    req(input$load_config$datapath)
    conf <- readRDS(input$load_config$datapath)
    
    # Restore file path and name
    file_path(conf$count_path)
    file_name(conf$file_name)
    meta_path(conf$meta_path)
    
    # Load count matrix
    if (!is.null(conf$count_path) && file.exists(conf$count_path)) {
      ext <- tools::file_ext(conf$count_path)
      count_data <- switch(ext,
                           csv = read.csv(conf$count_path, row.names = 1),
                           tsv = read.delim(conf$count_path, row.names = 1),
                           rds = readRDS(conf$count_path),
                           NULL)
      if (!is.null(count_data)) uploaded_data(count_data)
    }
    
    # Load metadata
    if (!is.null(conf$meta_path) && file.exists(conf$meta_path)) {
      ext <- tools::file_ext(conf$meta_path)
      meta_data <- switch(ext,
                          csv = read.csv(conf$meta_path, row.names = 1),
                          tsv = read.delim(conf$meta_path, row.names = 1),
                          rds = readRDS(conf$meta_path),
                          NULL)
      if (!is.null(meta_data)) metadata(meta_data)
    }
    
    showNotification("Data and metadata loaded from config", type = "message")
  })
  
  
  
  
  
  
  
  #### Upload Data
  
  shinyFileChoose(input, "user_file", roots = roots, filetypes = c("csv", "tsv", "rds"))
  
  # Fix path extraction: use [1,1] to grab actual string
  selected_path <- reactive({
    req(input$user_file)
    parseFilePaths(roots = roots, input$user_file)$datapath[1]
  })
  
  
  
  
  observeEvent(selected_path(), {
    req(selected_path())
    path <- selected_path()
    ext <- tools::file_ext(path)
    
    file_path(path)
    file_name(basename(path))
    
    if (!file.exists(path)) {
      showNotification("Count file not found. Check path or permissions.", type = "error")
    } else {
      count_data <- switch(ext,
                           csv = read.csv(path, row.names = 1),
                           tsv = read.delim(path, row.names = 1),
                           rds = {
                             dt <- readRDS(path)
                             dt %>% as.matrix() %>% as.data.frame()
                           },
                           NULL)
      
      if (!is.null(count_data)) {
        uploaded_data(count_data)
      
      }
    }
  })
  
  
  
  
  
  output$file_info <- renderPrint({
    req(selected_path())
    path <- selected_path()
    ext <- tools::file_ext(path)
    
    if (!file.exists(path)) {
      cat("File not found at:", path)
      return()
    }
    
    size_bytes <- file.info(path)$size
    size_mb <- size_bytes / 1024^2
    size_str <- if (size_mb > 1024) {
      paste0(round(size_mb / 1024, 2), " GB")
    } else {
      paste0(round(size_mb, 2), " MB")
    }
    
    df <- uploaded_data()
    cat("File path:", path, "\n")
    cat("File size:", size_str, "\n")
    cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
    cat("File type:", toupper(ext), "\n")
  })
  
  
  ## header file name 
  output$file_name_label <- renderText({
    if (is.null(file_name())) {
      return("Data not selected")
    }
    paste("Loaded:", file_name())
  })
  
  
  
  shinyFileChoose(input, "meta_file", roots = roots, filetypes = c("csv", "tsv", "rds"))
  
  
  

  
  
  metadata_path <- reactive({
    req(input$meta_file)
    parseFilePaths(roots = roots, input$meta_file)$datapath[1]
  })
  
  metadata <- reactive({
    req(config$meta_path)
    path <- config$meta_path
    ext <- tools::file_ext(path)
    
    meta_path(path)
    meta_file_name(basename(path))
    
    if (!file.exists(path)) {
      showNotification("Metadata file not found.", type = "error")
      return(NULL)
    }
    
    switch(ext,
           csv = read.csv(path, row.names = 1),
           tsv = read.delim(path, row.names = 1),
           rds = {
             dt <- readRDS(path)
             if (is.data.frame(dt)) return(dt)
             tryCatch(as.data.frame(as.matrix(dt)), error = function(e) NULL)
           },
           {
             showNotification("Unsupported metadata file type", type = "error")
             return(NULL)
           })
  })
  
  
  output$file_info <- renderPrint({
    req(selected_path())
    path <- selected_path()
    ext <- tools::file_ext(path)
    
    if (!file.exists(path)) {
      cat("File not found at:", path)
      return()
    }
    
    size_bytes <- file.info(path)$size
    size_mb <- size_bytes / 1024^2
    size_str <- if (size_mb > 1024) {
      paste0(round(size_mb / 1024, 2), " GB")
    } else {
      paste0(round(size_mb, 2), " MB")
    }
    
    df <- uploaded_data()
    cat("File path:", path, "\n")
    cat("File size:", size_str, "\n")
    cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
    cat("File type:", toupper(ext), "\n")
    
    # metadata info if available
    if (!is.null(input$meta_file)) {
      meta_path <- try(parseFilePaths(roots = roots, input$meta_file)$datapath[1], silent = TRUE)
      if (!inherits(meta_path, "try-error") && file.exists(meta_path)) {
        cat("Metadata file:", basename(meta_path), "\n")
      }
    }
  })
  
  
  
  ## gene tab server
  
  gene_annotation_result <- eventReactive(input$annotate_genes, {
    req(uploaded_data())
    
    gene_ids <- rownames(uploaded_data())
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = gene_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID")
    )
    
    # Join back to preserve order if needed
    annotated <- data.frame(ENSEMBL = gene_ids) %>%
      left_join(mapping, by = "ENSEMBL")
    
    annotated
  })
  
  output$gene_annotation_table <- DT::renderDataTable({
    req(gene_annotation_result())
    DT::datatable(gene_annotation_result(), options = list(pageLength = 10), rownames = FALSE)
  })
  
  
  
  
  observeEvent(input$convert_to_symbol, {
    req(uploaded_data())
    
    current_data <- uploaded_data()
    gene_ids <- rownames(current_data)
    
    mapping <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = gene_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    ) %>%
      distinct(ENSEMBL, SYMBOL) %>%
      filter(!is.na(SYMBOL)) %>%
      as.data.frame()
    
    mapping <- mapping[!duplicated(mapping$ENSEMBL), ]  # Remove duplicates manually
    
    common_ids <- intersect(rownames(current_data), mapping$ENSEMBL)
    
    new_data <- current_data[common_ids, , drop = FALSE]
    new_symbols <- mapping$SYMBOL[match(common_ids, mapping$ENSEMBL)]
    
    new_symbols <- make.unique(new_symbols)  # Ensure no duplicate rownames
    
    rownames(new_data) <- new_symbols
    
    uploaded_data(new_data)
    
    showNotification("Converted ENSEMBL IDs to Gene Symbols.", type = "message")
  })
  
  
  
  
  
  # gene norm
  log_normalize_seurat_style <- function(counts, scale_factor = 10000) {
    norm_counts <- t(t(counts) / colSums(counts)) * scale_factor
    log_norm <- log1p(norm_counts)
    return(log_norm)
  }
  
  get_hvg_vst_style <- function(counts, n_top = 2000, span = 0.3) {
    gene_means <- rowMeans(counts)
    gene_vars <- apply(counts, 1, var)
    names(gene_means) <- names(gene_vars) <- rownames(counts)
    
    keep <- gene_means > 0
    gene_means <- gene_means[keep]
    gene_vars <- gene_vars[keep]
    
    loess_fit <- loess(gene_vars ~ gene_means, span = span)
    expected_vars <- predict(loess_fit, newdata = gene_means)
    
    residual_var <- (gene_vars - expected_vars) / expected_vars
    residual_var[!is.finite(residual_var)] <- NA
    
    ranked_genes <- sort(residual_var, decreasing = TRUE, na.last = NA)
    top_genes <- names(ranked_genes)[1:min(n_top, length(ranked_genes))]
    
    return(list(top_genes = top_genes,
                mean = gene_means,
                sd = sqrt(gene_vars)))
  }
  
  
  
  ## hvg egenes
  hvg_result <- eventReactive(input$run_hvg, {
    req(uploaded_data())
    mat <- uploaded_data()
    
    # Filter genes with total counts >= 1
    mat <- mat[rowSums(mat) >= 1, , drop = FALSE]
    
    # Log-normalize
    log_norm <- log_normalize_seurat_style(mat)
    
    # HVG selection
    hvg_info <- get_hvg_vst_style(log_norm, n_top = input$top_n_hvg)
    
    # Optional overwrite
    if (input$overwrite_data) {
      updated <- uploaded_data()[hvg_info$top_genes, , drop = FALSE]
      uploaded_data(updated)
      
      
    } else if (input$overwrite_log_norm) {
      log_norm <- log_normalize_seurat_style(uploaded_data())
      uploaded_data(log_norm)
    }
    
    return(hvg_info)
  })
  
  
  
  output$saved_hvg_plot <- renderPlotly({
    req(hvg_plot())
    hvg_plot()
  })
  
  
  # plotly mean vs sd plot
  
  output$hvg_plot <- renderPlotly({
    req(hvg_result())
    gene_means <- hvg_result()$mean
    gene_sds <- hvg_result()$sd
    highlight <- rep("Other", length(gene_means))
    names(highlight) <- names(gene_means)
    highlight[hvg_result()$top_genes] <- "HVG"
    
    selected_row <- input$hvg_table_rows_selected
    selected_gene <- if (!is.null(selected_row) && length(selected_row) > 0) {
      hvg_result()$top_genes[selected_row]
    } else {
      NULL
    }
    
    color <- highlight
    if (!is.null(selected_gene)) {
      color[] <- "Other"
      color[selected_gene] <- "Selected"
    }
    
    df <- data.frame(
      Gene = names(gene_means),
      Mean = gene_means,
      SD = gene_sds,
      Group = factor(color, levels = c("Other", "HVG", "Selected"))
    )
    
    plot_obj <- plot_ly(
      data = df,
      x = ~Mean,
      y = ~SD,
      type = 'scatter',
      mode = 'markers',
      color = ~Group,
      colors = c("gray", "blue", "red"),
      text = ~paste0("Gene: ", Gene, "<br>Mean: ", round(Mean, 3), "<br>SD: ", round(SD, 3)),
      hoverinfo = 'text',
      marker = list(size = 5, opacity = 0.6)
    ) %>%
      layout(
        title = "HVG Selection: Mean vs SD",
        xaxis = list(title = "Mean (log-normalized)"),
        yaxis = list(title = "Standard Deviation")
      )
    
    hvg_plot(plot_obj)  # <-- Save it for reuse
    
    plot_obj
    
    
  })
  
  
  
  
  #hvg table
  output$hvg_table <- DT::renderDataTable({
    req(hvg_result())
    res <- hvg_result()
    top_genes <- res$top_genes
    gene_means <- res$mean[top_genes]
    gene_sds <- res$sd[top_genes]
    
    df <- data.frame(
      Gene = top_genes,
      Mean = round(gene_means, 4),
      SD = round(gene_sds, 4)
    )
    
    # allows to select a observation in the table
    #  input$hvg_table_rows_selected is automatically created by Shiny when you use DT::renderDataTable() with selection = "single" or "multiple"
    DT::datatable(df, selection = "single", options = list(pageLength = 10), rownames = FALSE)
  })
  
  
  # event with the selected observation (gene)
  observeEvent(input$hvg_table_rows_selected, {
    req(hvg_result())
    res <- hvg_result()
    selected_idx <- input$hvg_table_rows_selected
    if (!is.null(selected_idx)) {
      selected_hvg_gene(res$top_genes[selected_idx])
    }
  })
  
  
  
  
  
  #downloader
  
  output$download_hvg_list <- downloadHandler(
    filename = function() paste0("HVG_top_", input$top_n_hvg, ".csv"),
    content = function(file) {
      write.csv(data.frame(Gene = hvg_result()$top_genes), file, row.names = FALSE)
    }
  )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## UMAP server
  
  umap_plot_data <- eventReactive(input$plot_umap, {
    req(metadata(), input$umap_x, input$umap_y, input$umap_color)
    
    df <- metadata()
    
    x_vals <- df[[input$umap_x]]
    y_vals <- df[[input$umap_y]]
    color_vals <- df[[input$umap_color]]
    
    if (input$var_type == "Numerical") {
      color_vals <- suppressWarnings(as.numeric(color_vals))
    } else {
      color_vals <- as.factor(color_vals)
    }
    
    list(
      df = df,
      x_vals = x_vals,
      y_vals = y_vals,
      color_vals = color_vals
    )
  })
  
  
  output$umap_plot <- renderPlotly({
    req(umap_plot_data())
    vals <- umap_plot_data()
    
    plot_obj <- plot_ly(
      data = vals$df,
      x = ~vals$x_vals,
      y = ~vals$y_vals,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 6, opacity = 0.8),
      color = ~vals$color_vals,
      text = ~paste0(
        "ID: ", rownames(vals$df), "<br>",
        input$umap_x, ": ", vals$x_vals, "<br>",
        input$umap_y, ": ", vals$y_vals, "<br>",
        input$umap_color, ": ", vals$color_vals
      ),
      hoverinfo = 'text'
    ) %>% layout(title = "")
    
    hvg_plot(plot_obj)  # optionally save
    
    return(plot_obj)
  })
  
  
  
  output$umap_x_ui <- renderUI({
    req(metadata())
    selectInput("umap_x", "X Axis:", choices = colnames(metadata()), selected = colnames(metadata())[1])
  })
  
  output$umap_y_ui <- renderUI({
    req(metadata())
    selectInput("umap_y", "Y Axis:", choices = colnames(metadata()), selected = colnames(metadata())[2])
  })
  
  output$umap_color_ui <- renderUI({
    req(metadata())
    selectInput("umap_color", "Color by:", choices = colnames(metadata()))
  })
  
  
  
  
  ### Heatmap section
  
  # Populate gene and metadata options once files are loaded
  observeEvent(uploaded_data(), {
    req(uploaded_data())
    updatePickerInput(session, "genes_heatmap", choices = rownames(uploaded_data()))
  })
  
  observeEvent(metadata(), {
    req(metadata())
    updateSelectInput(session, "group_by_heatmap", choices = colnames(metadata()))
  })
  
  
  heatmap_matrix <- reactiveVal(NULL)
  
  
  
  
  # Generate Heatmap
  observeEvent(input$generate_heatmap, {
    output$heatmap_plot <- renderPlot({
      req(uploaded_data(), metadata(), input$genes_heatmap, input$group_by_heatmap)
      
      selected_genes <- input$genes_heatmap
      group_var <- input$group_by_heatmap
      
      
      
      
      log_norm <- uploaded_data()
      meta <- metadata()
      
      
      meta[[group_var]] <- as.factor(meta[[group_var]])
      
      desired_order <- unique(meta[[group_var]])
      
      custom_colors <- colorRamp2(c(-2, 0, 2), c("#008080", "#F5F5DC", "#FF7F50"))
      
      circle_fun <- function(x, y, w, h, v, color) {
        grid.rect(x = x, y = y, width = w, height = h, gp = gpar(fill = "white", col = NA))
        grid.circle(x = x, y = y, r = unit(v / 5000, "npc"), gp = gpar(fill = color, col = "black"))
      }
      
      log_norm_matrix <- log_norm[selected_genes, , drop = FALSE] %>%
        t() %>% as.data.frame() %>%
        mutate(group = meta[[group_var]][match(rownames(.), rownames(meta))])
      
      dt_matrix_percent <- log_norm[selected_genes, , drop = FALSE] %>%
        t() %>% as.data.frame() %>%
        mutate(group = meta[[group_var]][match(rownames(.), rownames(meta))]) %>%
        group_by(group) %>%
        summarize(across(all_of(selected_genes), ~ mean(.x > 0) * 100))
      
      log_norm_matrix <- log_norm_matrix %>%
        group_by(group) %>%
        summarize(across(all_of(selected_genes), mean, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(across(where(is.numeric), ~ ifelse(sd(.) == 0 | all(is.na(.)), 0, as.numeric(scale(.)))))
      
      mat_numeric <- t(as.matrix(log_norm_matrix[, -1]))
      colnames(mat_numeric) <- log_norm_matrix$group
      rownames(mat_numeric) <- selected_genes
      
      percent_matrix <- t(as.matrix(dt_matrix_percent[, -1]))
      colnames(percent_matrix) <- dt_matrix_percent$group
      rownames(percent_matrix) <- selected_genes
      
      mat_numeric <- mat_numeric[, desired_order, drop = FALSE]
      percent_matrix <- percent_matrix[, desired_order, drop = FALSE]
      
      heatmap_obj <- Heatmap(
        mat_numeric, name = "Scaled Expression", col = custom_colors,
        cell_fun = function(j, i, x, y, width, height, fill) {
          circle_fun(x, y, width, height, percent_matrix[i, j], fill)
        },
        show_row_names = TRUE, show_column_names = TRUE, cluster_columns = FALSE,
        column_names_gp = gpar(fontface = "bold", fontsize = 16, col = "black", rot = 45),
        row_names_gp = gpar(fontface = "bold", fontsize = 16, col = "black"),
        heatmap_legend_param = list(ticks_gp = gpar(fontface = "bold", fontsize = 16, col = "black")),
        column_names_rot = 45
      )
      
      heatmap_matrix(mat_numeric)
      
      draw(heatmap_obj)
      
      
      
      gene_expression(heatmap_obj)
      
      
    })
  })
  
  output$download_heatmap_matrix <- downloadHandler(
    filename = function() {
      paste0("Heatmap_Scaled_Matrix.xlsx")
    },
    content = function(file) {
      req(heatmap_matrix())
      write_xlsx(as.data.frame(heatmap_matrix()), path = file)
    }
  )
  
  
  
  
  
  
  ### Expression timeline part 
  
  
  
  observeEvent(uploaded_data(), {
    req(uploaded_data())
    updatePickerInput(session, "timeline_genes", choices = rownames(uploaded_data()))
  })
  
  observeEvent(metadata(), {
    req(metadata())
    updateSelectInput(session, "pseudo_var", choices = colnames(metadata()))
  })
  
  
  
  
  
  output$pseudotime_selector <- renderUI({
    req(metadata())
    selectInput("pseudotime_col", "Select Pseudotime Variable:", choices = colnames(metadata()))
  })
  
  output$timeline_heatmap <- renderPlot({
    req(uploaded_data(), metadata(), input$pseudo_var, input$timeline_genes)
    
    mat <- uploaded_data()
    meta <- metadata()
    
    pseudo_var <- input$pseudo_var
    bins <- input$timeline_bins
    genes <- input$timeline_genes
    
    # Ensure pseudotime column is numeric
    pseudotime_vals <- suppressWarnings(as.numeric(meta[[pseudo_var]]))
    if (anyNA(pseudotime_vals)) {
      showNotification("Selected pseudotime variable must be numeric.", type = "error")
      return(NULL)
    }
    
    # Prepare expression matrix
    matrix_heat <- t(mat[genes, , drop = FALSE]) %>% as.data.frame()
    matrix_heat$cell <- rownames(matrix_heat)
    
    # Add pseudotime
    meta <- meta %>% mutate(cell = rownames(meta))
    matrix_heat <- matrix_heat %>%
      left_join(meta %>% dplyr::select(cell, pseudotime = !!sym(pseudo_var)), by = "cell") %>%
      drop_na(pseudotime)
    
    # Bin pseudotime
    matrix_heat <- matrix_heat %>%
      arrange(pseudotime) %>%
      mutate(time_bin = ntile(pseudotime, bins))
    
    # Average expression per bin
    avg_expr <- matrix_heat %>%
      group_by(time_bin) %>%
      summarize(across(all_of(genes), mean, na.rm = TRUE)) %>%
      column_to_rownames("time_bin") %>%
      t() %>%
      as.data.frame()
    
    # Scale expression
    scaled_expr <- t(scale(t(avg_expr)))
    scaled_expr[scaled_expr > 2] <- 2
    scaled_expr[scaled_expr < -2] <- -2
    
    # Define color functions
    col_fun_expr <- colorRamp2(c(-2, 0, 2), scico(3, palette = "roma", direction = -1))
    col_fun_bins <- colorRamp2(1:bins, viridis::plasma(bins))
    
    bottom_anno <- HeatmapAnnotation("Pseudotime Bin" = 1:bins,
                                     col = list("Pseudotime Bin" = col_fun_bins),
                                     show_annotation_name = FALSE)
    
    draw(Heatmap(scaled_expr,
                 col = col_fun_expr,
                 name = "Scaled expression",
                 cluster_columns = FALSE,
                 cluster_rows = TRUE,
                 show_column_names = FALSE,
                 bottom_annotation = bottom_anno,
                 row_names_gp = gpar(fontsize = 10)))
    
    
    timeline_gene_expression(Heatmap(scaled_expr,
                                     col = col_fun_expr,
                                     name = "Scaled expression",
                                     cluster_columns = FALSE,
                                     cluster_rows = TRUE,
                                     show_column_names = FALSE,
                                     bottom_annotation = bottom_anno,
                                     row_names_gp = gpar(fontsize = 10)))
    
  })
  
  
  
  
  ##tissue
  
  observeEvent(metadata(), {
    updateSelectInput(session, "tissue_x", choices = colnames(metadata()))
    updateSelectInput(session, "tissue_y", choices = colnames(metadata()))
  })
  
  observeEvent(uploaded_data(), {
    updateSelectInput(session, "tissue_gene", choices = rownames(uploaded_data()))
  })
  
  
  
  output$tissue_plot <- renderPlotly({
    req(uploaded_data(), metadata(), input$tissue_x, input$tissue_y, input$tissue_gene)
    
    expr_matrix <- uploaded_data()
    meta <- metadata()
    gene <- input$tissue_gene
    
    # Get matching cell IDs
    cell_ids <- colnames(expr_matrix)
    meta_matched <- meta[cell_ids, , drop = FALSE]
    
    # Skip if gene not found or not enough data
    if (!(gene %in% rownames(expr_matrix))) return(NULL)
    
    # Expression for selected gene
    expr <- expr_matrix[gene, ]
    expr_log <- expr # disabled for now, it already log norma in the gene tab step!!
    expr_scaled <- as.numeric(scale(expr_log))
    
    df <- meta_matched %>%
      mutate(
        expr = expr_scaled,
        x = .data[[input$tissue_x]],
        y = .data[[input$tissue_y]]
      ) %>%
      drop_na(x, y, expr)
    
    tissue_plot_obj =  plot_ly(
      data = df,
      x = ~as.numeric(x),
      y = ~as.numeric(y),
      color = ~expr,
      type = "scatter",
      mode = "markers",
      colors = "viridis",
      marker = list(size = 6, opacity = 0.8),
      text = ~paste0("Expr: ", round(expr, 2)),
      hoverinfo = "text"
    ) %>% layout(title = paste("Expression of", gene))
    
    last_tissue_plot <<- tissue_plot_obj
    
    
    return(tissue_plot_obj)
  })
  
  
  observeEvent(input$save_tissue_to_dashboard, {
    req(last_tissue_plot)
    current <- tissue_dashboard_plots()
    tissue_dashboard_plots(c(current, list(last_tissue_plot)))
  })
  
  
  
  output$dashboard_tissue <- renderUI({
    plots <- tissue_dashboard_plots()
    if (length(plots) == 0) return(NULL)
    
    plot_output_list <- lapply(seq_along(plots), function(i) {
      plotname <- paste0("tissue_dashboard_", i)
      output[[plotname]] <- renderPlotly({ plots[[i]] })
      plotlyOutput(plotname)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  
  
  
  
  
  
  
  
  
  #dash
  
  
  output$saved_timeline_plot <- renderPlot({
    req(timeline_gene_expression())
    draw(timeline_gene_expression())
  })
  
  
  # interactive veiwer  save to umap dashboard
  observeEvent(input$save_umap_to_dashboard, {
    req(umap_plot_data())
    vals <- umap_plot_data()
    
    new_plot <- plot_ly(
      data = vals$df,
      x = ~vals$x_vals,
      y = ~vals$y_vals,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 6, opacity = 0.8),
      color = ~vals$color_vals,
      text = ~paste0(
        "ID: ", rownames(vals$df), "<br>",
        input$umap_x, ": ", vals$x_vals, "<br>",
        input$umap_y, ": ", vals$y_vals, "<br>",
        input$umap_color, ": ", vals$color_vals
      ),
      hoverinfo = 'text'
    ) %>% layout(title = paste("", input$umap_color))
    
    current_list <- umap_dashboard_plots()
    umap_dashboard_plots(c(current_list, list(new_plot)))
  })
  
  
  output$dashboard_umaps <- renderUI({
    plots <- umap_dashboard_plots()
    if (length(plots) == 0) return(h4("No plots saved yet."))
    
    plot_output_list <- lapply(seq_along(plots), function(i) {
      plotname <- paste0("umap_dashboard_", i)
      output[[plotname]] <- renderPlotly({ plots[[i]] })
      plotlyOutput(plotname)
    })
    
    do.call(tagList, plot_output_list)
  })
  
}




shinyApp(ui, server)

