# PPGLomics TCGA Explorer
# TCGA-PCPG Analysis Platform
# Developed by Hussam Alkaissi, MD MS
# With assistance from Claude AI

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)
library(DT)
library(RColorBrewer)
library(viridis)
library(grid)
library(ggrepel)
library(matrixTests)

options(warn = -1, shiny.maxRequestSize = 100*1024^2)

# ============================================================================
# LOAD DATA
# ============================================================================

tcga_expr <- read.csv("data/tcga_expression_log2.csv", row.names = 1, check.names = FALSE)
tcga_meta <- read.csv("data/tcga_metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)
colnames(tcga_meta) <- make.names(colnames(tcga_meta))

tcga_samples <- intersect(colnames(tcga_expr), tcga_meta$Sample.ID)
tcga_expr <- tcga_expr[, tcga_samples]
tcga_meta <- tcga_meta[match(tcga_samples, tcga_meta$Sample.ID), ]
rownames(tcga_meta) <- tcga_meta$Sample.ID
tcga_matrix <- as.matrix(tcga_expr)

# --- Gene lists ---
tcga_genes <- sort(rownames(tcga_expr))
tcga_genes_upper <- setNames(tcga_genes, toupper(tcga_genes))
tcga_gene_vars <- apply(tcga_matrix, 1, var, na.rm = TRUE)

# ============================================================================
# CONSTANTS & COLOR PALETTES
# ============================================================================

tcga_cluster_order <- c("Pseudohypoxia", "Kinase signaling", "Wnt-altered", "Cortical admixture")
tcga_genotype_order <- c("SDHx", "VHL", "EPAS1", "IDH1", "EGLN1", "NF1", "HRAS", "RET", 
                         "FGFR1", "NGFR", "BRAF", "SETD2", "MAML3 fusion", "CSDE1", 
                         "MAX", "TMEM127", "ARNT", "Sporadic/Unknown")
tcga_methylation_order <- c("low-methylated", "intermediate", "hyper-methylated")

tcga_genotype_colors <- c(
  "SDHx"="#E41A1C","VHL"="#377EB8","EPAS1"="#4DAF4A","IDH1"="#984EA3",
  "EGLN1"="#FF7F00","NF1"="#FFFF33","HRAS"="#A65628","RET"="#F781BF",
  "FGFR1"="#999999","NGFR"="#66C2A5","BRAF"="#FC8D62","SETD2"="#8DA0CB",
  "MAML3 fusion"="#E78AC3","CSDE1"="#A6D854","MAX"="#FFD92F",
  "TMEM127"="#E5C494","ARNT"="#B3B3B3","Sporadic/Unknown"="#7570B3")
tcga_cluster_colors <- c("Pseudohypoxia"="#E41A1C","Kinase signaling"="#377EB8",
                         "Wnt-altered"="#4DAF4A","Cortical admixture"="#984EA3")
tcga_metastatic_colors <- c("Yes"="#E41A1C","No"="#377EB8")
tcga_location_colors <- c("Adrenal (Pheo)"="#1B9E77","Extra-adrenal (PGL)"="#D95F02")
tcga_methylation_colors <- c("low-methylated"="#377EB8","intermediate"="#FFFF33","hyper-methylated"="#E41A1C")

gene_sets <- list(
  "Hypoxia/HIF"=c("VEGFA","SLC2A1","LDHA","PDK1","BNIP3","CA9","EPO","EGLN1","EGLN3","HK2","PGK1","ENO1"),
  "Chromaffin"=c("TH","DDC","DBH","PNMT","SLC6A2","SLC18A1","SLC18A2","CHGA","CHGB","NPY","PHOX2B","ASCL1"),
  "TCA Cycle"=c("CS","ACO2","IDH2","IDH3A","OGDH","SUCLA2","SDHA","SDHB","SDHC","SDHD","FH","MDH2"),
  "Imprinted"=c("IGF2","H19","DLK1","MEG3","MEST","PEG3","PEG10","SNRPN","NDN","PLAGL1","CDKN1C"),
  "Kinase"=c("RET","NF1","HRAS","KRAS","NRAS","BRAF","RAF1","MAP2K1","MAPK1","MAPK3","FGFR1","EGFR"),
  "Wnt"=c("WNT1","WNT3A","WNT5A","CTNNB1","APC","AXIN1","AXIN2","GSK3B","LEF1","TCF7","MYC","CCND1")
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

find_genes <- function(input_genes, gene_list, gene_upper) {
  input_genes <- trimws(input_genes)
  input_genes <- input_genes[input_genes != ""]
  found <- c()
  for (g in input_genes) {
    if (g %in% gene_list) found <- c(found, g)
    else if (toupper(g) %in% names(gene_upper)) found <- c(found, gene_upper[[toupper(g)]])
  }
  unique(found)
}

remove_outliers <- function(x, coef = 1.5) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  x[x < q[1] - coef * iqr | x > q[2] + coef * iqr] <- NA
  x
}

# ============================================================================
# UI
# ============================================================================

ui <- navbarPage(
  title = tags$span(
    tags$strong("PPGLomics TCGA-PCPG", style = "color: white;"),
    tags$span(" | ", style = "color: #888; margin: 0 5px;"),
    tags$span("Alkaissi Lab", style = "color: #E8654B; font-weight: 500;")
  ),
  id = "main_nav",
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  tags$head(tags$style(HTML("
    .btn-secondary, .btn-default {
      background-color: #2C3E50 !important;
      border-color: #2C3E50 !important;
      color: white !important;
    }
    .btn-secondary:hover, .btn-default:hover {
      background-color: #1a252f !important;
      border-color: #1a252f !important;
    }
    .download-btn {
      background-color: #2C3E50 !important;
      border-color: #2C3E50 !important;
      color: white !important;
      margin-bottom: 5px;
    }
    .dataset-info {
      background-color: #3498db;
      color: white;
      padding: 8px 15px;
      border-radius: 5px;
      font-weight: bold;
      margin-bottom: 15px;
    }
  "))),
  
  # ============================================================================
  # SINGLE GENE TAB
  # ============================================================================
  tabPanel("Single Gene", value = "single_gene",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "TCGA-PCPG (n=160)"),
        selectizeInput("sg_gene", "Gene:", choices = NULL, 
                       options = list(placeholder = "Type gene name...", maxOptions = 100)),
        selectInput("sg_group", "Group by:", 
                    c("Primary Genotype"="Primary_Genotype", 
                      "mRNA Cluster"="mRNA.Subtype.Clusters",
                      "Metastatic Status"="Metastatic", 
                      "Tumor Location"="Location",
                      "Methylation Cluster"="Methylation.Cluster")),
        hr(),
        h5("Display Options"),
        checkboxInput("sg_points", "Show points", TRUE),
        checkboxInput("sg_violin", "Violin plot", FALSE),
        checkboxInput("sg_remove_outliers", "Remove outliers", FALSE),
        hr(),
        downloadButton("sg_dl", "Download Plot", class = "download-btn"),
        downloadButton("sg_dl_csv", "Download Data (CSV)", class = "download-btn")
      ),
      mainPanel(width = 9, 
        plotOutput("sg_plot", height = "500px"), 
        verbatimTextOutput("sg_stats")
      )
    )
  ),
  
  # ============================================================================
  # HEATMAP TAB
  # ============================================================================
  tabPanel("Heatmap", value = "heatmap",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "TCGA-PCPG (n=160)"),
        selectInput("hm_preset", "Preset Gene Set:", c("Select..."="", names(gene_sets))),
        textAreaInput("hm_genes", "Genes (one per line):", rows = 8, placeholder = "PNMT\nth\nDbh"),
        actionButton("hm_go", "Generate Heatmap", class = "btn-primary"),
        hr(),
        selectInput("hm_annot", "Annotations:", multiple = TRUE,
                    c("Genotype"="Primary_Genotype", "Cluster"="mRNA.Subtype.Clusters",
                      "Metastatic"="Metastatic", "Location"="Location"),
                    selected = c("Primary_Genotype", "mRNA.Subtype.Clusters")),
        checkboxInput("hm_scale", "Z-score scale", TRUE),
        checkboxInput("hm_clust_row", "Cluster rows", TRUE),
        checkboxInput("hm_clust_col", "Cluster columns", TRUE),
        downloadButton("hm_dl", "Download Heatmap", class = "download-btn")
      ),
      mainPanel(width = 9, plotOutput("hm_plot", height = "650px"))
    )
  ),
  
  # ============================================================================
  # CORRELATION TAB
  # ============================================================================
  tabPanel("Correlation", value = "correlation",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "TCGA-PCPG (n=160)"),
        h4("Two-Gene Correlation"),
        selectizeInput("cor_g1", "Gene 1:", choices = NULL, 
                       options = list(placeholder = "Type gene name...")),
        selectizeInput("cor_g2", "Gene 2:", choices = NULL,
                       options = list(placeholder = "Type gene name...")),
        selectInput("cor_color", "Color by:", 
                    c("None"="none", "Genotype"="Primary_Genotype", "Cluster"="mRNA.Subtype.Clusters",
                      "Metastatic"="Metastatic", "Location"="Location")),
        selectInput("cor_method", "Method:", c("Spearman"="spearman", "Pearson"="pearson")),
        checkboxInput("cor_show_line", "Show regression line", TRUE),
        hr(),
        h4("Find Similar Genes"),
        selectizeInput("sim_gene", "Reference gene:", choices = NULL,
                       options = list(placeholder = "Type gene name...")),
        numericInput("sim_n", "Top N:", 25, 5, 100),
        actionButton("sim_go", "Find Similar", class = "btn-primary"),
        downloadButton("cor_dl", "Download Plot", class = "download-btn"),
        downloadButton("cor_dl_csv", "Download Similar Genes (CSV)", class = "download-btn")
      ),
      mainPanel(width = 9,
        plotOutput("cor_plot", height = "350px"),
        verbatimTextOutput("cor_stats"),
        fluidRow(
          column(6, h4("Positively Correlated"), DTOutput("sim_pos")),
          column(6, h4("Negatively Correlated"), DTOutput("sim_neg"))
        )
      )
    )
  ),
  
  # ============================================================================
  # VOLCANO PLOT TAB
  # ============================================================================
  tabPanel("Volcano Plot", value = "volcano",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "TCGA-PCPG (n=160)"),
        selectInput("vol_var", "Compare by:", 
                    c("mRNA Cluster"="mRNA.Subtype.Clusters", "Genotype"="Primary_Genotype",
                      "Metastatic Status"="Metastatic", "Tumor Location"="Location",
                      "Methylation Pattern"="Methylation.Cluster")),
        uiOutput("vol_g1_ui"), 
        uiOutput("vol_g2_ui"),
        hr(),
        numericInput("vol_fc", "log2FC threshold:", 1, 0, 5, 0.25),
        numericInput("vol_p", "-log10(p) threshold:", 2, 0, 10, 0.5),
        textAreaInput("vol_hl", "Highlight genes:", rows = 2, placeholder = "PNMT\nEPAS1"),
        checkboxInput("vol_labels", "Label top genes", TRUE),
        conditionalPanel("input.vol_labels",
          checkboxInput("vol_label_fc", "By fold change", TRUE),
          checkboxInput("vol_label_pval", "By p-value", TRUE),
          numericInput("vol_nlabel", "Labels per direction:", 5, 0, 15)
        ),
        hr(),
        actionButton("vol_go", "Run Analysis", class = "btn-primary"),
        downloadButton("vol_dl_plot", "Download Plot", class = "download-btn"),
        downloadButton("vol_dl_data", "Download Results (CSV)", class = "download-btn")
      ),
      mainPanel(width = 9,
        h4("Volcano Plot"),
        plotOutput("vol_plot", height = "550px"),
        DTOutput("vol_table")
      )
    )
  ),
  
  # ============================================================================
  # SURVIVAL TAB
  # ============================================================================
  tabPanel("Survival", value = "survival",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "TCGA-PCPG (n=160)"),
        selectInput("surv_ep", "Endpoint:", c("Overall Survival"="os", "Event-Free (Aggressive)"="efs")),
        selectInput("surv_by", "Stratify by:", 
                    c("Genotype"="Primary_Genotype", "Cluster"="mRNA.Subtype.Clusters",
                      "Metastatic"="Metastatic", "Location"="Location", "Gene Expression"="gene")),
        conditionalPanel("input.surv_by == 'gene'",
          selectizeInput("surv_gene", "Gene:", choices = NULL,
                         options = list(placeholder = "Type gene name...")),
          selectInput("surv_split", "Split:", c("Median"="median", "Tertile"="tertile"))
        ),
        hr(),
        actionButton("surv_go", "Run Analysis", class = "btn-primary"),
        downloadButton("surv_dl", "Download Plot", class = "download-btn")
      ),
      mainPanel(width = 9, 
        plotOutput("surv_plot", height = "500px"), 
        verbatimTextOutput("surv_stats")
      )
    )
  ),
  
  # ============================================================================
  # ONCOPRINT TAB
  # ============================================================================
  tabPanel("OncoPrint", value = "oncoprint",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "TCGA-PCPG (n=160)"),
        h4("Driver Alterations"),
        p("Red = Somatic | Blue = Germline | Green = Fusion"),
        hr(),
        downloadButton("onco_dl", "Download OncoPrint", class = "download-btn")
      ),
      mainPanel(width = 9, 
        plotOutput("onco_plot", height = "600px")
      )
    )
  ),
  
  # ============================================================================
  # ABOUT TAB
  # ============================================================================
  tabPanel("About", value = "about",
    fluidRow(
      column(8, offset = 2, style = "padding-top: 30px;",
        h2("PPGLomics TCGA Explorer", style = "color: #2C3E50;"),
        p("Transcriptomic analysis platform for the TCGA-PCPG dataset.", 
          style = "font-size: 18px; color: #666;"),
        hr(),
        
        h4("Dataset & Reference", style = "color: #2C3E50;"),
        p(strong("TCGA-PCPG:"), " n=160 samples, multiple genotypes"),
        p("Fishbein L, et al. Comprehensive Molecular Characterization of Pheochromocytoma and Paraganglioma. ",
          em("Cancer Cell."), " 2017;31(2):181-193. ",
          tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/28162975", "PMID: 28162975", target = "_blank"),
          style = "font-size: 13px; color: #666; margin-left: 15px;"),
        
        hr(),
        h4("Developer", style = "color: #2C3E50;"),
        p(tags$strong("Hussam Alkaissi, MD MS")),
        p("Clinician-scientist investigating and managing patients with pheochromocytoma and paraganglioma, 
           with research interests in pseudohypoxia, hypoxia signaling, and Krebs cycle defects."),
        p(em("Developed with assistance from Claude AI")),
        
        hr(),
        p(em("PPGLomics TCGA v1.0 | Alkaissi Lab"), style = "color: #999; text-align: center;")
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # Update gene selectors
  observe({
    updateSelectizeInput(session, "sg_gene", choices = tcga_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "cor_g1", choices = tcga_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "cor_g2", choices = tcga_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "sim_gene", choices = tcga_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "surv_gene", choices = tcga_genes, selected = "", server = TRUE)
  })

  # ==========================================================================
  # SINGLE GENE ANALYSIS
  # ==========================================================================
  
  sg_data <- reactive({
    req(input$sg_gene, input$sg_group)
    req(input$sg_gene %in% tcga_genes, input$sg_gene != "")
    
    gv <- input$sg_group
    df <- data.frame(Sample = colnames(tcga_matrix), Expr = as.numeric(tcga_matrix[input$sg_gene, ]))
    df <- merge(df, tcga_meta, by.x = "Sample", by.y = "Sample.ID")
    
    if (isTRUE(input$sg_remove_outliers)) {
      df$Expr <- remove_outliers(df$Expr)
      df <- df[!is.na(df$Expr), ]
    }
    
    df <- df[!is.na(df[[gv]]) & df[[gv]] != "",]
    
    # Filter TCGA specific
    if (gv == "Metastatic") df <- df[df[[gv]] %in% c("Yes", "No"), ]
    if (gv == "Location") df <- df[df[[gv]] %in% c("Adrenal (Pheo)", "Extra-adrenal (PGL)"), ]
    
    order_map <- list(
      "Metastatic" = c("No", "Yes"),
      "Location" = c("Adrenal (Pheo)", "Extra-adrenal (PGL)"),
      "mRNA.Subtype.Clusters" = tcga_cluster_order,
      "Methylation.Cluster" = tcga_methylation_order,
      "Primary_Genotype" = tcga_genotype_order
    )
    
    if (gv %in% names(order_map)) {
      lvls <- intersect(order_map[[gv]], unique(df[[gv]]))
      df[[gv]] <- factor(df[[gv]], levels = lvls)
    } else {
      df[[gv]] <- factor(df[[gv]])
    }
    df
  })
  
  sg_plot_obj <- reactive({
    df <- sg_data()
    req(nrow(df) > 0)
    gv <- input$sg_group
    
    ns <- df %>% count(.data[[gv]])
    df <- left_join(df, ns, by = gv)
    df$lbl <- paste0(df[[gv]], "\n(n=", df$n, ")")
    df$lbl <- factor(df$lbl, levels = unique(df$lbl[order(as.numeric(df[[gv]]))]))
    
    cols <- switch(gv, "Primary_Genotype"=tcga_genotype_colors, "mRNA.Subtype.Clusters"=tcga_cluster_colors,
                   "Metastatic"=tcga_metastatic_colors, "Location"=tcga_location_colors, 
                   "Methylation.Cluster"=tcga_methylation_colors, NULL)
    
    p <- ggplot(df, aes(lbl, Expr, fill = .data[[gv]])) + 
      theme_classic(base_size = 14) +
      labs(x = NULL, y = paste0(input$sg_gene, " Expression (log2 TPM+1)"), title = input$sg_gene) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5))
    
    if (input$sg_violin) p <- p + geom_violin(alpha = 0.6, scale = "width")
    p <- p + geom_boxplot(width = ifelse(input$sg_violin, 0.2, 0.5), outlier.shape = NA, alpha = 0.8)
    
    if (input$sg_points) p <- p + geom_jitter(width = 0.15, size = 1.5, alpha = 0.5)
    
    if (!is.null(cols)) p <- p + scale_fill_manual(values = cols)
    p
  })
  
  output$sg_plot <- renderPlot({ req(input$sg_gene != ""); sg_plot_obj() })
  
  output$sg_stats <- renderPrint({
    req(input$sg_gene != "")
    df <- sg_data(); req(nrow(df) > 0)
    gv <- input$sg_group
    kw <- kruskal.test(df$Expr ~ df[[gv]])
    cat("Kruskal-Wallis p =", format(kw$p.value, digits = 3), "\n\n")
    print(df %>% group_by(.data[[gv]]) %>% 
            summarise(n=n(), mean=round(mean(Expr),2), median=round(median(Expr),2), .groups="drop") %>% 
            as.data.frame())
  })
  
  output$sg_dl <- downloadHandler(
    filename = function() paste0(input$sg_gene, "_TCGA_expression.pdf"),
    content = function(file) { pdf(file, width = 10, height = 7); print(sg_plot_obj()); dev.off() }
  )
  
  output$sg_dl_csv <- downloadHandler(
    filename = function() paste0(input$sg_gene, "_TCGA_data.csv"),
    content = function(file) {
      df <- sg_data()
      gv <- input$sg_group
      out <- df[, c("Sample", "Expr", gv)]
      colnames(out) <- c("Sample_ID", paste0(input$sg_gene, "_Expression"), "Group")
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # HEATMAP
  # ==========================================================================
  
  observeEvent(input$hm_preset, {
    if (input$hm_preset %in% names(gene_sets)) {
      genes_available <- intersect(gene_sets[[input$hm_preset]], tcga_genes)
      updateTextAreaInput(session, "hm_genes", value = paste(genes_available, collapse = "\n"))
    }
  })
  
  hm_g <- eventReactive(input$hm_go, {
    find_genes(strsplit(input$hm_genes, "\n")[[1]], tcga_genes, tcga_genes_upper)
  }, ignoreNULL = FALSE)
  
  hm_plot_obj <- reactive({
    genes <- hm_g(); req(length(genes) >= 2)
    mat <- tcga_matrix[genes, , drop = FALSE]
    if (input$hm_scale) mat <- t(scale(t(mat)))
    
    ann_cols <- list(Primary_Genotype = tcga_genotype_colors, mRNA.Subtype.Clusters = tcga_cluster_colors,
                     Metastatic = tcga_metastatic_colors, Location = tcga_location_colors)
    
    ann_df <- tcga_meta[colnames(mat), input$hm_annot, drop = FALSE]
    ann_cols_filtered <- ann_cols[intersect(names(ann_cols), input$hm_annot)]
    
    ha <- HeatmapAnnotation(df = ann_df, col = ann_cols_filtered, na_col = "grey90")
    rng <- max(abs(mat), na.rm = TRUE)
    col_fun <- colorRamp2(c(-rng, 0, rng), c("#2166AC", "white", "#B2182B"))
    
    Heatmap(mat, name = ifelse(input$hm_scale, "Z-score", "Expression"), col = col_fun, 
            top_annotation = ha, cluster_rows = input$hm_clust_row, cluster_columns = input$hm_clust_col,
            show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),
            column_title = paste0(length(genes), " genes × ", ncol(mat), " samples"))
  })
  
  output$hm_plot <- renderPlot({ hm_plot_obj() })
  
  output$hm_dl <- downloadHandler(
    filename = function() "heatmap_TCGA.pdf",
    content = function(file) { pdf(file, 14, 10); draw(hm_plot_obj()); dev.off() }
  )

  # ==========================================================================
  # CORRELATION
  # ==========================================================================
  
  cor_plot_obj <- reactive({
    req(input$cor_g1 %in% tcga_genes, input$cor_g2 %in% tcga_genes)
    req(input$cor_g1 != "", input$cor_g2 != "")
    
    g1 <- as.numeric(tcga_matrix[input$cor_g1,]); g2 <- as.numeric(tcga_matrix[input$cor_g2,])
    df <- data.frame(x = g1, y = g2, Sample = colnames(tcga_matrix))
    df <- merge(df, tcga_meta, by.x = "Sample", by.y = "Sample.ID")
    ct <- cor.test(g1, g2, method = input$cor_method)
    
    p <- ggplot(df, aes(x, y))
    
    if (!is.null(input$cor_color) && input$cor_color != "none") {
      p <- p + geom_point(aes_string(color = input$cor_color), size = 2.5, alpha = 0.7)
      cols <- switch(input$cor_color, "Primary_Genotype"=tcga_genotype_colors, 
                     "mRNA.Subtype.Clusters"=tcga_cluster_colors, "Metastatic"=tcga_metastatic_colors,
                     "Location"=tcga_location_colors, NULL)
      if (!is.null(cols)) p <- p + scale_color_manual(values = cols)
    } else {
      p <- p + geom_point(size = 2.5, alpha = 0.7, color = "#377EB8")
    }
    
    if (isTRUE(input$cor_show_line)) {
      p <- p + geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE)
    }
    
    p + labs(x = paste0(input$cor_g1, " (log2 TPM+1)"), y = paste0(input$cor_g2, " (log2 TPM+1)"),
             title = paste(input$cor_g1, "vs", input$cor_g2),
             subtitle = sprintf("%s r = %.3f, p = %.2e", tools::toTitleCase(input$cor_method), ct$estimate, ct$p.value)) + 
      theme_classic(base_size = 14)
  })
  
  output$cor_plot <- renderPlot({ req(input$cor_g1 != "", input$cor_g2 != ""); cor_plot_obj() })
  
  output$cor_stats <- renderPrint({
    req(input$cor_g1 %in% tcga_genes, input$cor_g2 %in% tcga_genes)
    req(input$cor_g1 != "", input$cor_g2 != "")
    ct <- cor.test(as.numeric(tcga_matrix[input$cor_g1,]), as.numeric(tcga_matrix[input$cor_g2,]), method = input$cor_method)
    cat(sprintf("%s r = %.4f, p = %.2e, n = %d\n", tools::toTitleCase(input$cor_method), ct$estimate, ct$p.value, ncol(tcga_matrix)))
  })
  
  sim_rv <- reactiveValues(pos = NULL, neg = NULL, query = NULL)
  
  observeEvent(input$sim_go, {
    req(input$sim_gene %in% tcga_genes, input$sim_gene != "")
    withProgress(message = "Computing correlations...", {
      target <- tcga_matrix[input$sim_gene, ]
      valid <- setdiff(names(tcga_gene_vars[tcga_gene_vars > 0]), input$sim_gene)
      cors <- cor(target, t(tcga_matrix[valid, ]), method = "spearman", use = "pairwise.complete.obs")[1, ]
      cors <- cors[!is.na(cors)]
      n <- min(input$sim_n, length(cors))
      pos <- sort(cors, decreasing = TRUE)[1:n]
      neg <- sort(cors)[1:n]
      sim_rv$pos <- data.frame(Gene = names(pos), r = round(unname(pos), 4))
      sim_rv$neg <- data.frame(Gene = names(neg), r = round(unname(neg), 4))
      sim_rv$query <- input$sim_gene
    })
  })
  
  output$sim_pos <- renderDT({ req(sim_rv$pos); datatable(sim_rv$pos, selection = "single", options = list(pageLength = 8, dom = "tp"), rownames = FALSE) })
  output$sim_neg <- renderDT({ req(sim_rv$neg); datatable(sim_rv$neg, selection = "single", options = list(pageLength = 8, dom = "tp"), rownames = FALSE) })
  
  observeEvent(input$sim_pos_rows_selected, {
    sel <- input$sim_pos_rows_selected
    if (length(sel) > 0 && !is.null(sim_rv$pos)) {
      updateSelectizeInput(session, "cor_g1", selected = sim_rv$query)
      updateSelectizeInput(session, "cor_g2", selected = sim_rv$pos$Gene[sel])
    }
  })
  
  observeEvent(input$sim_neg_rows_selected, {
    sel <- input$sim_neg_rows_selected
    if (length(sel) > 0 && !is.null(sim_rv$neg)) {
      updateSelectizeInput(session, "cor_g1", selected = sim_rv$query)
      updateSelectizeInput(session, "cor_g2", selected = sim_rv$neg$Gene[sel])
    }
  })
  
  output$cor_dl <- downloadHandler(
    filename = function() paste0(input$cor_g1, "_vs_", input$cor_g2, "_TCGA.pdf"),
    content = function(file) { pdf(file, 8, 6); print(cor_plot_obj()); dev.off() }
  )
  
  output$cor_dl_csv <- downloadHandler(
    filename = function() paste0("similar_genes_", input$sim_gene, "_TCGA.csv"),
    content = function(file) {
      req(sim_rv$pos, sim_rv$neg)
      pos <- sim_rv$pos; pos$Direction <- "Positive"
      neg <- sim_rv$neg; neg$Direction <- "Negative"
      out <- rbind(pos, neg)
      write.csv(out, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # VOLCANO PLOT
  # ==========================================================================
  
  vol_choices <- reactive({
    req(input$vol_var)
    switch(input$vol_var, 
           "mRNA.Subtype.Clusters"=tcga_cluster_order,
           "Primary_Genotype"=intersect(tcga_genotype_order, unique(tcga_meta$Primary_Genotype)),
           "Metastatic"=c("Yes", "No"),
           "Location"=c("Adrenal (Pheo)","Extra-adrenal (PGL)"),
           "Methylation.Cluster"=tcga_methylation_order)
  })
  
  output$vol_g1_ui <- renderUI({
    ch <- vol_choices()
    selectInput("vol_g1", "Group 1 (numerator):", ch, multiple = TRUE, selected = ch[1])
  })
  
  output$vol_g2_ui <- renderUI({
    ch <- vol_choices()
    remaining <- setdiff(ch, input$vol_g1)
    if (length(remaining) == 0) remaining <- ch
    selectInput("vol_g2", "Group 2 (denominator):", remaining, multiple = TRUE, selected = remaining[1])
  })
  
  vol_rv <- reactiveValues(data = NULL, l1 = NULL, l2 = NULL, n1 = NULL, n2 = NULL)
  
  observeEvent(input$vol_go, {
    req(input$vol_g1, input$vol_g2)
    col <- input$vol_var
    
    s1 <- tcga_meta$Sample.ID[tcga_meta[[col]] %in% input$vol_g1]
    s2 <- tcga_meta$Sample.ID[tcga_meta[[col]] %in% input$vol_g2]
    s1 <- intersect(s1, colnames(tcga_matrix)); s2 <- intersect(s2, colnames(tcga_matrix))
    req(length(s1) >= 3, length(s2) >= 3)
    
    withProgress(message = "Running differential expression...", {
      m1 <- tcga_matrix[, s1, drop = FALSE]; m2 <- tcga_matrix[, s2, drop = FALSE]
      res <- row_t_welch(m1, m2)
      res$gene <- rownames(res)
      res$log2FC <- res$mean.x - res$mean.y
      res$neglog10p <- -log10(res$pvalue)
      res$padj <- p.adjust(res$pvalue, method = "BH")
      vol_rv$data <- res
      vol_rv$l1 <- paste(input$vol_g1, collapse = "+")
      vol_rv$l2 <- paste(input$vol_g2, collapse = "+")
      vol_rv$n1 <- length(s1); vol_rv$n2 <- length(s2)
    })
  })
  
  vol_plot_obj <- reactive({
    req(vol_rv$data)
    df <- vol_rv$data
    df$sig <- ifelse(abs(df$log2FC) >= input$vol_fc & df$neglog10p >= input$vol_p, "Significant", "NS")
    df$direction <- ifelse(df$log2FC > 0, "Upregulated", "Downregulated")
    df$color_group <- ifelse(df$sig == "NS", "NS", df$direction)
    
    hl_genes <- find_genes(strsplit(input$vol_hl, "\n")[[1]], tcga_genes, tcga_genes_upper)
    df$highlight <- df$gene %in% hl_genes
    
    n_per_dir <- input$vol_nlabel
    top_genes <- data.frame()
    
    if (isTRUE(input$vol_label_fc)) {
      top_fc_up <- df %>% filter(log2FC > 0) %>% arrange(desc(abs(log2FC))) %>% head(n_per_dir)
      top_fc_down <- df %>% filter(log2FC < 0) %>% arrange(desc(abs(log2FC))) %>% head(n_per_dir)
      top_genes <- rbind(top_genes, top_fc_up, top_fc_down)
    }
    
    if (isTRUE(input$vol_label_pval)) {
      top_p_up <- df %>% filter(log2FC > 0) %>% arrange(desc(neglog10p)) %>% head(n_per_dir)
      top_p_down <- df %>% filter(log2FC < 0) %>% arrange(desc(neglog10p)) %>% head(n_per_dir)
      top_genes <- rbind(top_genes, top_p_up, top_p_down)
    }
    
    top_genes <- top_genes[!duplicated(top_genes$gene), ]
    
    color_map <- c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "NS" = "grey70")
    
    p <- ggplot(df, aes(log2FC, neglog10p)) +
      geom_point(aes(color = color_group), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = color_map, name = "Direction",
                         breaks = c("Upregulated", "Downregulated", "NS"),
                         labels = c("Upregulated", "Downregulated", "Not Significant")) +
      geom_vline(xintercept = c(-input$vol_fc, input$vol_fc), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = input$vol_p, linetype = "dashed", color = "grey50") +
      theme_classic(base_size = 14) +
      labs(x = "log2 Fold Change", y = "-log10(p-value)",
           title = paste0(vol_rv$l1, " (n=", vol_rv$n1, ") vs ", vol_rv$l2, " (n=", vol_rv$n2, ")")) +
      theme(legend.position = "right")
    
    if (input$vol_labels && nrow(top_genes) > 0) {
      p <- p + geom_text_repel(data = top_genes, aes(label = gene), 
                               size = 3, max.overlaps = 30, segment.color = "grey50")
    }
    
    if (any(df$highlight)) {
      hl_df <- df[df$highlight, ]
      p <- p + geom_point(data = hl_df, color = "black", size = 3, shape = 1, stroke = 1.5) +
        geom_text_repel(data = hl_df, aes(label = gene), size = 3.5, fontface = "bold",
                        color = "black", max.overlaps = 30)
    }
    p
  })
  
  output$vol_plot <- renderPlot({ req(vol_rv$data); vol_plot_obj() })
  
  output$vol_table <- renderDT({
    req(vol_rv$data)
    df <- vol_rv$data %>%
      filter(abs(log2FC) >= input$vol_fc & neglog10p >= input$vol_p) %>%
      select(gene, log2FC, pvalue, padj) %>% arrange(pvalue) %>%
      mutate(log2FC = round(log2FC, 3), pvalue = signif(pvalue, 3), padj = signif(padj, 3))
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$vol_dl_plot <- downloadHandler(
    filename = function() "volcano_TCGA.pdf",
    content = function(file) { pdf(file, 10, 8); print(vol_plot_obj()); dev.off() }
  )
  
  output$vol_dl_data <- downloadHandler(
    filename = function() "DEG_results_TCGA.csv",
    content = function(file) {
      req(vol_rv$data)
      write.csv(vol_rv$data %>% select(gene, log2FC, pvalue, padj) %>% arrange(pvalue), file, row.names = FALSE)
    }
  )

  # ==========================================================================
  # SURVIVAL
  # ==========================================================================
  
  surv_rv <- reactiveValues(data = NULL, title = NULL)
  
  observeEvent(input$surv_go, {
    df <- tcga_meta
    
    if (input$surv_ep == "os") {
      df$time <- as.numeric(df$Max.Days.of.Follow.up.from.Initial.Surgery) / 30.44
      df$status <- ifelse(df$Vital.Status == "Dead", 1, 0)
    } else {
      df$time <- as.numeric(df$Event.Free.Days..Aggressive.disease.) / 30.44
      df$status <- ifelse(df$Clinically.Aggressive.and.or.Metastatic %in% c("Yes", "yes"), 1, 0)
    }
    
    df <- df[!is.na(df$time) & df$time > 0 & !is.na(df$status), ]
    
    if (input$surv_by == "gene") {
      req(input$surv_gene %in% tcga_genes, input$surv_gene != "")
      df <- df[df$Sample.ID %in% colnames(tcga_matrix), ]
      ev <- as.numeric(tcga_matrix[input$surv_gene, df$Sample.ID])
      if (input$surv_split == "median") {
        df$group <- ifelse(ev >= median(ev, na.rm = TRUE), "High", "Low")
      } else {
        qs <- quantile(ev, c(1/3, 2/3), na.rm = TRUE)
        df$group <- cut(ev, c(-Inf, qs, Inf), labels = c("Low", "Mid", "High"))
      }
      surv_rv$title <- paste0("Survival by ", input$surv_gene, " Expression")
    } else {
      df$group <- df[[input$surv_by]]
      surv_rv$title <- paste0("Survival by ", input$surv_by)
    }
    
    df <- df[!is.na(df$group) & df$group != "" & df$group != "Unknown", ]
    req(nrow(df) >= 10)
    surv_rv$data <- df
  })
  
  surv_plot_obj <- reactive({
    req(surv_rv$data)
    df <- surv_rv$data
    fit <- survfit(Surv(time, status) ~ group, data = df)
    lr <- survdiff(Surv(time, status) ~ group, data = df)
    pv <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
    
    ggsurvplot(fit, data = df, conf.int = FALSE, risk.table = TRUE, 
               xlab = "Months", ggtheme = theme_classic(),
               title = surv_rv$title, subtitle = sprintf("Log-rank p = %.2e", pv))
  })
  
  output$surv_plot <- renderPlot({ req(surv_rv$data); print(surv_plot_obj()) })
  
  output$surv_stats <- renderPrint({
    req(surv_rv$data)
    df <- surv_rv$data
    lr <- survdiff(Surv(time, status) ~ group, data = df)
    pv <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
    cat("Log-rank p =", format(pv, digits = 3), "\n\nGroups:\n")
    print(table(df$group))
  })
  
  output$surv_dl <- downloadHandler(
    filename = function() "survival_TCGA.pdf",
    content = function(file) {
      req(surv_rv$data)
      p <- surv_plot_obj()
      ggsave(file, plot = p$plot, width = 10, height = 8)
    }
  )

  # ==========================================================================
  # ONCOPRINT
  # ==========================================================================
  
  onco_plot_obj <- reactive({
    genes <- c("SDHx","VHL","EPAS1","IDH1","EGLN1","NF1","HRAS","RET","FGFR1","NGFR","BRAF","SETD2","MAML3 fusion","CSDE1","MAX","TMEM127","ARNT")
    mat <- matrix("", nrow = length(genes), ncol = nrow(tcga_meta))
    rownames(mat) <- genes; colnames(mat) <- tcga_meta$Sample.ID
    
    has_mut <- function(x) !is.na(x) & x != 0 & x != "0" & x != ""
    
    for (gene in genes) {
      if (gene == "SDHx") {
        for (i in seq_len(nrow(tcga_meta)))
          if (has_mut(tcga_meta$SDHB.Germline.Mutation[i]) || has_mut(tcga_meta$SDHD.Germline.Mutation[i]))
            mat["SDHx", tcga_meta$Sample.ID[i]] <- "Germline"
      } else if (gene == "MAML3 fusion") {
        for (i in seq_len(nrow(tcga_meta)))
          if (isTRUE(tcga_meta$UBTF.MAML3_fusion[i] == 1) || isTRUE(tcga_meta$X.MAML3_fusion[i] == 1))
            mat["MAML3 fusion", tcga_meta$Sample.ID[i]] <- "Fusion"
      } else if (gene == "NGFR") {
        for (i in seq_len(nrow(tcga_meta)))
          if (isTRUE(tcga_meta$X.NGFR_fusion[i] == 1)) mat["NGFR", tcga_meta$Sample.ID[i]] <- "Fusion"
      } else if (gene %in% c("MAX","TMEM127","EGLN1")) {
        col <- paste0(gene, ".Germline.Mutation")
        if (col %in% names(tcga_meta))
          for (i in seq_len(nrow(tcga_meta)))
            if (has_mut(tcga_meta[[col]][i])) mat[gene, tcga_meta$Sample.ID[i]] <- "Germline"
      } else {
        scol <- paste0(gene, ".Somatic.Mutation"); gcol <- paste0(gene, ".Germline.Mutation")
        for (i in seq_len(nrow(tcga_meta))) {
          if (scol %in% names(tcga_meta) && has_mut(tcga_meta[[scol]][i])) mat[gene, tcga_meta$Sample.ID[i]] <- "Somatic"
          if (gcol %in% names(tcga_meta) && has_mut(tcga_meta[[gcol]][i])) mat[gene, tcga_meta$Sample.ID[i]] <- "Germline"
        }
      }
    }
    
    mat <- mat[order(rowSums(mat != ""), decreasing = TRUE), ]
    mat <- mat[, order(colSums(mat != ""), decreasing = TRUE)]
    mat <- mat[, colSums(mat != "") > 0, drop = FALSE]
    
    freq <- round(100 * rowSums(mat != "") / ncol(mat), 1)
    alter_fun <- list(background = alter_graphic("rect", fill = "#EEEEEE"),
                      Somatic = alter_graphic("rect", fill = "#E41A1C"),
                      Germline = alter_graphic("rect", fill = "#377EB8"),
                      Fusion = alter_graphic("rect", fill = "#4DAF4A"))
    
    oncoPrint(mat, alter_fun = alter_fun, col = c(Somatic = "#E41A1C", Germline = "#377EB8", Fusion = "#4DAF4A"),
              show_column_names = FALSE, alter_fun_is_vectorized = FALSE,
              right_annotation = rowAnnotation(Freq = anno_barplot(freq, width = unit(2.5, "cm"))),
              column_title = paste0("TCGA-PCPG Driver Alterations (", ncol(mat), " samples)"))
  })
  
  output$onco_plot <- renderPlot({ onco_plot_obj() })
  
  output$onco_dl <- downloadHandler(
    filename = function() "oncoprint_TCGA.pdf",
    content = function(file) { pdf(file, 14, 8); draw(onco_plot_obj()); dev.off() }
  )
}

shinyApp(ui, server)
