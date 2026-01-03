# PPGLomics A5 Explorer
# A5 SDHB Consortium Analysis Platform
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

a5_expr <- read.csv("data/a5_expression_log2.csv", row.names = 1, check.names = FALSE)
a5_meta <- read.csv("data/a5_metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)

a5_samples <- intersect(colnames(a5_expr), a5_meta$SAMPLE_ID)
a5_expr <- a5_expr[, a5_samples]
a5_meta <- a5_meta[match(a5_samples, a5_meta$SAMPLE_ID), ]
rownames(a5_meta) <- a5_meta$SAMPLE_ID
a5_matrix <- as.matrix(a5_expr)

# --- Fix A5 labels ---
a5_meta$TUMOR_SIZE <- gsub("^gt5cm$", ">5cm", a5_meta$TUMOR_SIZE)
a5_meta$TUMOR_SIZE <- gsub("^le5cm$", "<=5cm", a5_meta$TUMOR_SIZE)
a5_meta$TUMOR_SIZE_DISPLAY <- gsub("<=5cm", "≤5cm", a5_meta$TUMOR_SIZE)
a5_meta$KI67_GROUP <- gsub("^Ki67_low$", "Ki67 <1%", a5_meta$KI67_GROUP)
a5_meta$KI67_GROUP <- gsub("^Ki67_1_5$", "Ki67 1-5%", a5_meta$KI67_GROUP)
a5_meta$KI67_GROUP <- gsub("^Ki67_5_10$", "Ki67 5-10%", a5_meta$KI67_GROUP)
a5_meta$KI67_GROUP <- gsub("^Ki67_high$", "Ki67 >10%", a5_meta$KI67_GROUP)
a5_meta$METASTATIC <- ifelse(a5_meta$LESION_TYPE == "Metastasis", "Metastatic", "Non-Metastatic")

# --- Gene lists ---
a5_genes <- sort(rownames(a5_expr))
a5_genes_upper <- setNames(a5_genes, toupper(a5_genes))
a5_gene_vars <- apply(a5_matrix, 1, var, na.rm = TRUE)

# ============================================================================
# CONSTANTS & COLOR PALETTES
# ============================================================================

a5_lesion_order <- c("Primary_No_Mets", "Primary_With_Mets", "Metastasis")
a5_location_order <- c("Adrenal", "HNPGL", "T_PGL", "AP_PGL", "Metastatic")
a5_biochem_order <- c("Non_Secreting", "NE", "NE_DA", "DA")
a5_tert_atrx_order <- c("WT", "TERT_altered", "ATRX_altered")
a5_mutation_order <- c("Missense", "Nonsense", "Frameshift", "Splice", "Deletion")
a5_grantham_order <- c("Conservative", "Moderately_Conservative", "Moderately_Radical", "Radical")
a5_ki67_order <- c("Ki67 <1%", "Ki67 1-5%", "Ki67 5-10%", "Ki67 >10%")
a5_tumor_size_order <- c("<=5cm", ">5cm")

a5_lesion_colors <- c("Primary_No_Mets"="#377EB8", "Primary_With_Mets"="#FF7F00", "Metastasis"="#E41A1C")
a5_location_colors <- c("Adrenal"="#1B9E77", "HNPGL"="#D95F02", "T_PGL"="#7570B3", 
                        "AP_PGL"="#E6AB02", "Metastatic"="#E41A1C")
a5_biochem_colors <- c("Non_Secreting"="#999999", "NE"="#E41A1C", "NE_DA"="#FF7F00", "DA"="#377EB8")
a5_tert_atrx_colors <- c("WT"="#377EB8", "TERT_altered"="#E41A1C", "ATRX_altered"="#4DAF4A")
a5_mutation_colors <- c("Missense"="#377EB8", "Nonsense"="#E41A1C", "Frameshift"="#4DAF4A", 
                        "Splice"="#984EA3", "Deletion"="#FF7F00")
a5_grantham_colors <- c("Conservative"="#377EB8", "Moderately_Conservative"="#4DAF4A",
                        "Moderately_Radical"="#FF7F00", "Radical"="#E41A1C")
a5_ki67_colors <- c("Ki67 <1%"="#377EB8", "Ki67 1-5%"="#4DAF4A", "Ki67 5-10%"="#FF7F00", "Ki67 >10%"="#E41A1C")
a5_chromothripsis_colors <- c("No"="#377EB8", "Yes"="#E41A1C")
a5_sex_colors <- c("Male"="#377EB8", "Female"="#E377C2")
a5_binary_colors <- c("Yes"="#E41A1C", "No"="#377EB8")
a5_tumor_size_colors <- c("<=5cm"="#377EB8", ">5cm"="#E41A1C")
a5_metastatic_colors <- c("Non-Metastatic"="#377EB8", "Metastatic"="#E41A1C")
a5_prior_treatment_colors <- c("Yes"="#E41A1C", "No"="#377EB8")

shape_values_hollow <- c("WT"=16, "TERT_altered"=2, "ATRX_altered"=0,
                         "Ki67 <1%"=16, "Ki67 1-5%"=2, "Ki67 5-10%"=0, "Ki67 >10%"=5,
                         "No"=16, "Yes"=2,
                         "Non-Metastatic"=16, "Metastatic"=2)

gene_sets <- list(
  "Hypoxia/HIF"=c("VEGFA","SLC2A1","LDHA","PDK1","BNIP3","CA9","EPO","EGLN1","EGLN3","HK2","PGK1","ENO1"),
  "Chromaffin"=c("TH","DDC","DBH","PNMT","SLC6A2","SLC18A1","SLC18A2","CHGA","CHGB","NPY","PHOX2B","ASCL1"),
  "TCA Cycle"=c("CS","ACO2","IDH2","IDH3A","OGDH","SUCLA2","SDHA","SDHB","SDHC","SDHD","FH","MDH2"),
  "Imprinted"=c("IGF2","H19","DLK1","MEG3","MEST","PEG3","PEG10","SNRPN","NDN","PLAGL1","CDKN1C"),
  "Kinase"=c("RET","NF1","HRAS","KRAS","NRAS","BRAF","RAF1","MAP2K1","MAPK1","MAPK3","FGFR1","EGFR"),
  "Wnt"=c("WNT1","WNT3A","WNT5A","CTNNB1","APC","AXIN1","AXIN2","GSK3B","LEF1","TCF7","MYC","CCND1")
)

a5_paired_info <- list(
  "E143" = list(primary = c("E143-P1"), metastasis = c("E143-M1", "E143-M2"), variant = "p.Arg90*"),
  "E146" = list(primary = c("E146-P1"), metastasis = c("E146-M1"), variant = "p.Arg46*"),
  "E158" = list(primary = c("E158-P1"), metastasis = c("E158-M1"), variant = "p.Val140Phe"),
  "E159" = list(primary = c("E159-P2"), metastasis = c("E159-M1", "E159-M2"), variant = "p.Val140Phe"),
  "E225" = list(primary = c("E225-P1"), metastasis = c("E225-M1"), variant = "p.Asp64Thrfs*13")
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
    tags$strong("PPGLomics A5", style = "color: white;"),
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
      background-color: #e74c3c;
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
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        selectizeInput("sg_gene", "Gene:", choices = NULL, 
                       options = list(placeholder = "Type gene name...", maxOptions = 100)),
        selectInput("sg_group", "Group by:", 
                    c("Lesion Type"="LESION_TYPE", "Location"="LOCATION_GROUP",
                      "TERT/ATRX Status"="TERT_ATRX_STATUS", "Mutation Type"="MUTATION_TYPE",
                      "Grantham Class"="GRANTHAM_CLASS", "Biochemical"="BIOCHEM_PHENOTYPE",
                      "Ki67 Group"="KI67_GROUP", "Tumor Size"="TUMOR_SIZE",
                      "Chromothripsis"="CHROMOTHRIPSIS", "Hypertension"="HTN", "Sex"="SEX",
                      "Prior DNA Damage"="PRIOR_TREATMENT")),
        hr(),
        h5("Display Options"),
        checkboxInput("sg_points", "Show points", TRUE),
        checkboxInput("sg_violin", "Violin plot", FALSE),
        checkboxInput("sg_remove_outliers", "Remove outliers", FALSE),
        hr(),
        h5("Overlay"),
        selectInput("sg_overlay", "Overlay by:", 
                    c("None"="none", "TERT/ATRX Status"="TERT_ATRX_STATUS",
                      "Ki67 Group"="KI67_GROUP", "Chromothripsis"="CHROMOTHRIPSIS",
                      "Metastatic"="METASTATIC")),
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
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        selectInput("hm_preset", "Preset Gene Set:", c("Select..."="", names(gene_sets))),
        textAreaInput("hm_genes", "Genes (one per line):", rows = 8, placeholder = "PNMT\nth\nDbh"),
        actionButton("hm_go", "Generate Heatmap", class = "btn-primary"),
        hr(),
        selectInput("hm_annot", "Annotations:", multiple = TRUE,
                    c("Lesion Type"="LESION_TYPE", "Location"="LOCATION_GROUP",
                      "TERT/ATRX"="TERT_ATRX_STATUS", "Mutation Type"="MUTATION_TYPE",
                      "Biochemical"="BIOCHEM_PHENOTYPE", "Ki67"="KI67_GROUP",
                      "Chromothripsis"="CHROMOTHRIPSIS", "Hypertension"="HTN",
                      "Prior DNA Damage"="PRIOR_TREATMENT"),
                    selected = c("LESION_TYPE", "TERT_ATRX_STATUS")),
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
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        h4("Two-Gene Correlation"),
        selectizeInput("cor_g1", "Gene 1:", choices = NULL, 
                       options = list(placeholder = "Type gene name...")),
        selectizeInput("cor_g2", "Gene 2:", choices = NULL,
                       options = list(placeholder = "Type gene name...")),
        selectInput("cor_color", "Color by:", 
                    c("None"="none", "Lesion Type"="LESION_TYPE", "Location"="LOCATION_GROUP",
                      "TERT/ATRX"="TERT_ATRX_STATUS", "Biochemical"="BIOCHEM_PHENOTYPE",
                      "Hypertension"="HTN", "Tumor Size"="TUMOR_SIZE")),
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
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        selectInput("vol_var", "Compare by:", 
                    c("Paired: Metastasis vs Primary"="PAIRED",
                      "Lesion Type"="LESION_TYPE", "Location"="LOCATION_GROUP",
                      "TERT/ATRX Status"="TERT_ATRX_STATUS", "Mutation Type"="MUTATION_TYPE",
                      "Grantham Class"="GRANTHAM_CLASS", "Biochemical"="BIOCHEM_PHENOTYPE",
                      "Ki67 Group"="KI67_GROUP", "Chromothripsis"="CHROMOTHRIPSIS",
                      "Hypertension"="HTN", "Sex"="SEX", "Prior DNA Damage"="PRIOR_TREATMENT")),
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
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        selectInput("surv_ep", "Endpoint:", c("Overall Survival"="os", "Metastasis-Free Survival"="mfs")),
        selectInput("surv_by", "Stratify by:", 
                    c("Lesion Type"="LESION_TYPE", "Location"="LOCATION_GROUP",
                      "TERT/ATRX Status"="TERT_ATRX_STATUS", "Mutation Type"="MUTATION_TYPE",
                      "Grantham Class"="GRANTHAM_CLASS", "Ki67 Group"="KI67_GROUP",
                      "Chromothripsis"="CHROMOTHRIPSIS", "Tumor Size"="TUMOR_SIZE",
                      "Hypertension"="HTN", "Gene Expression"="gene")),
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
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        h4("Genomic Alterations"),
        h5("Alteration Rows"),
        checkboxInput("onco_tert", "TERT promoter", TRUE),
        checkboxInput("onco_atrx", "ATRX loss", TRUE),
        checkboxInput("onco_chromo", "Chromothripsis", TRUE),
        hr(),
        h5("Annotation Tracks"),
        checkboxInput("onco_lesion", "Lesion Type", TRUE),
        checkboxInput("onco_location", "Location", FALSE),
        checkboxInput("onco_mutation", "Mutation Type", TRUE),
        checkboxInput("onco_htn", "Hypertension", FALSE),
        checkboxInput("onco_sex", "Sex", FALSE),
        checkboxInput("onco_size", "Tumor Size", FALSE),
        checkboxInput("onco_ki67", "Ki67", FALSE),
        checkboxInput("onco_grantham", "Grantham Class", FALSE),
        checkboxInput("onco_prior", "Prior DNA Damage", FALSE),
        hr(),
        downloadButton("onco_dl", "Download OncoPrint", class = "download-btn")
      ),
      mainPanel(width = 9, 
        plotOutput("onco_plot", height = "600px"), 
        DTOutput("onco_table")
      )
    )
  ),
  
  # ============================================================================
  # PAIRED ANALYSIS TAB
  # ============================================================================
  tabPanel("Paired Analysis", value = "paired",
    sidebarLayout(
      sidebarPanel(width = 3,
        div(class = "dataset-info", "A5 SDHB Consortium (n=91)"),
        h4("Primary vs Metastasis"),
        p("5 patients with paired samples", style = "color: #666; font-size: 12px;"),
        selectizeInput("paired_gene", "Gene:", choices = NULL,
                       options = list(placeholder = "Type gene name...")),
        checkboxInput("paired_show_all", "Show all paired patients", TRUE),
        conditionalPanel("!input.paired_show_all",
          selectInput("paired_patient", "Select Patient:", 
                      choices = c("E143", "E146", "E158", "E159", "E225"))
        ),
        hr(),
        actionButton("paired_go", "Generate Expression Plot", class = "btn-primary"),
        downloadButton("paired_dl", "Download Plot", class = "download-btn"),
        downloadButton("paired_dl_csv", "Download Data (CSV)", class = "download-btn")
      ),
      mainPanel(width = 9,
        tabsetPanel(id = "paired_tabs",
          tabPanel("Expression Plot",
            plotOutput("paired_plot", height = "450px"),
            verbatimTextOutput("paired_stats")
          ),
          tabPanel("Patient Details",
            DTOutput("paired_info_table")
          )
        )
      )
    )
  ),
  
  # ============================================================================
  # ABOUT TAB
  # ============================================================================
  tabPanel("About", value = "about",
    fluidRow(
      column(8, offset = 2, style = "padding-top: 30px;",
        h2("PPGLomics A5 Explorer", style = "color: #2C3E50;"),
        p("Transcriptomic analysis platform for the A5 SDHB Consortium dataset.", 
          style = "font-size: 18px; color: #666;"),
        hr(),
        
        h4("Dataset & Reference", style = "color: #2C3E50;"),
        p(strong("A5 Consortium:"), " n=91 samples, SDHB germline mutations only"),
        p("Flynn A, et al. Multi-omic analysis of SDHB-deficient pheochromocytomas and paragangliomas identifies metastasis and treatment-related molecular profiles. ",
          em("Nat Commun."), " 2025;16(1):2632. ",
          tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/40097403", "PMID: 40097403", target = "_blank"),
          style = "font-size: 13px; color: #666; margin-left: 15px;"),
        
        hr(),
        h4("Developer", style = "color: #2C3E50;"),
        p(tags$strong("Hussam Alkaissi, MD MS")),
        p("Clinician-scientist investigating and managing patients with pheochromocytoma and paraganglioma, 
           with research interests in pseudohypoxia, hypoxia signaling, and Krebs cycle defects."),
        p(em("Developed with assistance from Claude AI")),
        
        hr(),
        p(em("PPGLomics A5 v1.0 | Alkaissi Lab"), style = "color: #999; text-align: center;")
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
    updateSelectizeInput(session, "sg_gene", choices = a5_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "cor_g1", choices = a5_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "cor_g2", choices = a5_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "sim_gene", choices = a5_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "surv_gene", choices = a5_genes, selected = "", server = TRUE)
    updateSelectizeInput(session, "paired_gene", choices = a5_genes, selected = "", server = TRUE)
  })

  # ==========================================================================
  # SINGLE GENE ANALYSIS
  # ==========================================================================
  
  sg_data <- reactive({
    req(input$sg_gene, input$sg_group)
    req(input$sg_gene %in% a5_genes, input$sg_gene != "")
    
    gv <- input$sg_group
    df <- data.frame(Sample = colnames(a5_matrix), Expr = as.numeric(a5_matrix[input$sg_gene, ]))
    df <- merge(df, a5_meta, by.x = "Sample", by.y = "SAMPLE_ID")
    
    if (isTRUE(input$sg_remove_outliers)) {
      df$Expr <- remove_outliers(df$Expr)
      df <- df[!is.na(df$Expr), ]
    }
    
    df <- df[!is.na(df[[gv]]) & df[[gv]] != "",]
    
    order_map <- list(
      "LESION_TYPE" = a5_lesion_order, "LOCATION_GROUP" = a5_location_order,
      "BIOCHEM_PHENOTYPE" = a5_biochem_order, "TERT_ATRX_STATUS" = a5_tert_atrx_order,
      "MUTATION_TYPE" = a5_mutation_order, "GRANTHAM_CLASS" = a5_grantham_order,
      "KI67_GROUP" = a5_ki67_order, "TUMOR_SIZE" = a5_tumor_size_order
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
    
    cols <- switch(gv, "LESION_TYPE"=a5_lesion_colors, "LOCATION_GROUP"=a5_location_colors,
                   "BIOCHEM_PHENOTYPE"=a5_biochem_colors, "TERT_ATRX_STATUS"=a5_tert_atrx_colors,
                   "MUTATION_TYPE"=a5_mutation_colors, "GRANTHAM_CLASS"=a5_grantham_colors,
                   "KI67_GROUP"=a5_ki67_colors, "CHROMOTHRIPSIS"=a5_chromothripsis_colors,
                   "SEX"=a5_sex_colors, "HTN"=a5_binary_colors, "TUMOR_SIZE"=a5_tumor_size_colors,
                   "PRIOR_TREATMENT"=a5_prior_treatment_colors, NULL)
    
    p <- ggplot(df, aes(lbl, Expr, fill = .data[[gv]])) + 
      theme_classic(base_size = 14) +
      labs(x = NULL, y = paste0(input$sg_gene, " Expression (log2 CPM)"), title = input$sg_gene) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5))
    
    if (input$sg_violin) p <- p + geom_violin(alpha = 0.6, scale = "width")
    p <- p + geom_boxplot(width = ifelse(input$sg_violin, 0.2, 0.5), outlier.shape = NA, alpha = 0.8)
    
    if (input$sg_points) {
      if (!is.null(input$sg_overlay) && input$sg_overlay != "none") {
        overlay_var <- input$sg_overlay
        df_overlay <- df[!is.na(df[[overlay_var]]) & df[[overlay_var]] != "", ]
        p <- p + geom_jitter(data = df_overlay, 
                             aes_string(shape = overlay_var), 
                             width = 0.15, size = 3, alpha = 0.8, stroke = 1.5) +
          scale_shape_manual(values = shape_values_hollow, name = overlay_var) +
          guides(fill = "none") +
          theme(legend.position = "right")
      } else {
        p <- p + geom_jitter(width = 0.15, size = 1.5, alpha = 0.5)
      }
    }
    
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
    filename = function() paste0(input$sg_gene, "_A5_expression.pdf"),
    content = function(file) { pdf(file, width = 10, height = 7); print(sg_plot_obj()); dev.off() }
  )
  
  output$sg_dl_csv <- downloadHandler(
    filename = function() paste0(input$sg_gene, "_A5_data.csv"),
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
      genes_available <- intersect(gene_sets[[input$hm_preset]], a5_genes)
      updateTextAreaInput(session, "hm_genes", value = paste(genes_available, collapse = "\n"))
    }
  })
  
  hm_g <- eventReactive(input$hm_go, {
    find_genes(strsplit(input$hm_genes, "\n")[[1]], a5_genes, a5_genes_upper)
  }, ignoreNULL = FALSE)
  
  hm_plot_obj <- reactive({
    genes <- hm_g(); req(length(genes) >= 2)
    mat <- a5_matrix[genes, , drop = FALSE]
    if (input$hm_scale) mat <- t(scale(t(mat)))
    
    ann_cols <- list(LESION_TYPE = a5_lesion_colors, LOCATION_GROUP = a5_location_colors,
                     TERT_ATRX_STATUS = a5_tert_atrx_colors, MUTATION_TYPE = a5_mutation_colors,
                     BIOCHEM_PHENOTYPE = a5_biochem_colors, KI67_GROUP = a5_ki67_colors,
                     CHROMOTHRIPSIS = a5_chromothripsis_colors, HTN = a5_binary_colors,
                     PRIOR_TREATMENT = a5_prior_treatment_colors)
    
    ann_df <- a5_meta[colnames(mat), input$hm_annot, drop = FALSE]
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
    filename = function() "heatmap_A5.pdf",
    content = function(file) { pdf(file, 14, 10); draw(hm_plot_obj()); dev.off() }
  )

  # ==========================================================================
  # CORRELATION
  # ==========================================================================
  
  cor_plot_obj <- reactive({
    req(input$cor_g1 %in% a5_genes, input$cor_g2 %in% a5_genes)
    req(input$cor_g1 != "", input$cor_g2 != "")
    
    g1 <- as.numeric(a5_matrix[input$cor_g1,]); g2 <- as.numeric(a5_matrix[input$cor_g2,])
    df <- data.frame(x = g1, y = g2, Sample = colnames(a5_matrix))
    df <- merge(df, a5_meta, by.x = "Sample", by.y = "SAMPLE_ID")
    ct <- cor.test(g1, g2, method = input$cor_method)
    
    p <- ggplot(df, aes(x, y))
    
    if (!is.null(input$cor_color) && input$cor_color != "none") {
      p <- p + geom_point(aes_string(color = input$cor_color), size = 2.5, alpha = 0.7)
      cols <- switch(input$cor_color, "LESION_TYPE"=a5_lesion_colors, "LOCATION_GROUP"=a5_location_colors,
                     "TERT_ATRX_STATUS"=a5_tert_atrx_colors, "BIOCHEM_PHENOTYPE"=a5_biochem_colors,
                     "HTN"=a5_binary_colors, "TUMOR_SIZE"=a5_tumor_size_colors, NULL)
      if (!is.null(cols)) p <- p + scale_color_manual(values = cols)
    } else {
      p <- p + geom_point(size = 2.5, alpha = 0.7, color = "#377EB8")
    }
    
    if (isTRUE(input$cor_show_line)) {
      p <- p + geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE)
    }
    
    p + labs(x = paste0(input$cor_g1, " (log2 CPM)"), y = paste0(input$cor_g2, " (log2 CPM)"),
             title = paste(input$cor_g1, "vs", input$cor_g2),
             subtitle = sprintf("%s r = %.3f, p = %.2e", tools::toTitleCase(input$cor_method), ct$estimate, ct$p.value)) + 
      theme_classic(base_size = 14)
  })
  
  output$cor_plot <- renderPlot({ req(input$cor_g1 != "", input$cor_g2 != ""); cor_plot_obj() })
  
  output$cor_stats <- renderPrint({
    req(input$cor_g1 %in% a5_genes, input$cor_g2 %in% a5_genes)
    req(input$cor_g1 != "", input$cor_g2 != "")
    ct <- cor.test(as.numeric(a5_matrix[input$cor_g1,]), as.numeric(a5_matrix[input$cor_g2,]), method = input$cor_method)
    cat(sprintf("%s r = %.4f, p = %.2e, n = %d\n", tools::toTitleCase(input$cor_method), ct$estimate, ct$p.value, ncol(a5_matrix)))
  })
  
  sim_rv <- reactiveValues(pos = NULL, neg = NULL, query = NULL)
  
  observeEvent(input$sim_go, {
    req(input$sim_gene %in% a5_genes, input$sim_gene != "")
    withProgress(message = "Computing correlations...", {
      target <- a5_matrix[input$sim_gene, ]
      valid <- setdiff(names(a5_gene_vars[a5_gene_vars > 0]), input$sim_gene)
      cors <- cor(target, t(a5_matrix[valid, ]), method = "spearman", use = "pairwise.complete.obs")[1, ]
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
    filename = function() paste0(input$cor_g1, "_vs_", input$cor_g2, "_A5.pdf"),
    content = function(file) { pdf(file, 8, 6); print(cor_plot_obj()); dev.off() }
  )
  
  output$cor_dl_csv <- downloadHandler(
    filename = function() paste0("similar_genes_", input$sim_gene, "_A5.csv"),
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
    if (input$vol_var == "PAIRED") return(NULL)
    switch(input$vol_var,
           "LESION_TYPE"=a5_lesion_order,
           "LOCATION_GROUP"=intersect(a5_location_order, unique(a5_meta$LOCATION_GROUP)),
           "TERT_ATRX_STATUS"=a5_tert_atrx_order,
           "MUTATION_TYPE"=intersect(a5_mutation_order, unique(a5_meta$MUTATION_TYPE)),
           "GRANTHAM_CLASS"=intersect(a5_grantham_order, unique(a5_meta$GRANTHAM_CLASS)),
           "BIOCHEM_PHENOTYPE"=intersect(a5_biochem_order, unique(a5_meta$BIOCHEM_PHENOTYPE)),
           "KI67_GROUP"=intersect(a5_ki67_order, unique(a5_meta$KI67_GROUP)),
           "CHROMOTHRIPSIS"=c("Yes", "No"),
           "HTN"=c("Yes", "No"),
           "SEX"=c("Male", "Female"),
           "PRIOR_TREATMENT"=c("Yes", "No"),
           unique(na.omit(a5_meta[[input$vol_var]])))
  })
  
  output$vol_g1_ui <- renderUI({
    req(input$vol_var)
    if (input$vol_var == "PAIRED") {
      helpText("Compares Metastasis vs Primary from 5 patients with paired samples")
    } else {
      ch <- vol_choices()
      selectInput("vol_g1", "Group 1 (numerator):", ch, multiple = TRUE, selected = ch[1])
    }
  })
  
  output$vol_g2_ui <- renderUI({
    req(input$vol_var)
    if (input$vol_var == "PAIRED") return(NULL)
    ch <- vol_choices()
    remaining <- setdiff(ch, input$vol_g1)
    if (length(remaining) == 0) remaining <- ch
    selectInput("vol_g2", "Group 2 (denominator):", remaining, multiple = TRUE, selected = remaining[1])
  })
  
  vol_rv <- reactiveValues(data = NULL, l1 = NULL, l2 = NULL, n1 = NULL, n2 = NULL, mode = NULL)
  
  observeEvent(input$vol_go, {
    if (!is.null(input$vol_var) && input$vol_var == "PAIRED") {
      primary_samples <- unlist(lapply(a5_paired_info, function(x) x$primary))
      met_samples <- unlist(lapply(a5_paired_info, function(x) x$metastasis))
      primary_samples <- intersect(primary_samples, colnames(a5_matrix))
      met_samples <- intersect(met_samples, colnames(a5_matrix))
      
      req(length(primary_samples) >= 2, length(met_samples) >= 2)
      
      withProgress(message = "Running paired differential expression...", {
        m1 <- a5_matrix[, met_samples, drop = FALSE]
        m2 <- a5_matrix[, primary_samples, drop = FALSE]
        res <- row_t_welch(m1, m2)
        res$gene <- rownames(res)
        res$log2FC <- res$mean.x - res$mean.y
        res$neglog10p <- -log10(res$pvalue)
        res$padj <- p.adjust(res$pvalue, method = "BH")
        vol_rv$data <- res
        vol_rv$l1 <- "Metastasis"
        vol_rv$l2 <- "Primary"
        vol_rv$n1 <- length(met_samples)
        vol_rv$n2 <- length(primary_samples)
        vol_rv$mode <- "paired"
      })
    } else {
      req(input$vol_g1, input$vol_g2)
      col <- input$vol_var
      
      s1 <- a5_meta$SAMPLE_ID[a5_meta[[col]] %in% input$vol_g1]
      s2 <- a5_meta$SAMPLE_ID[a5_meta[[col]] %in% input$vol_g2]
      s1 <- intersect(s1, colnames(a5_matrix)); s2 <- intersect(s2, colnames(a5_matrix))
      req(length(s1) >= 3, length(s2) >= 3)
      
      withProgress(message = "Running differential expression...", {
        m1 <- a5_matrix[, s1, drop = FALSE]; m2 <- a5_matrix[, s2, drop = FALSE]
        res <- row_t_welch(m1, m2)
        res$gene <- rownames(res)
        res$log2FC <- res$mean.x - res$mean.y
        res$neglog10p <- -log10(res$pvalue)
        res$padj <- p.adjust(res$pvalue, method = "BH")
        vol_rv$data <- res
        vol_rv$l1 <- paste(input$vol_g1, collapse = "+")
        vol_rv$l2 <- paste(input$vol_g2, collapse = "+")
        vol_rv$n1 <- length(s1); vol_rv$n2 <- length(s2)
        vol_rv$mode <- "standard"
      })
    }
  })
  
  vol_plot_obj <- reactive({
    req(vol_rv$data)
    df <- vol_rv$data
    df$sig <- ifelse(abs(df$log2FC) >= input$vol_fc & df$neglog10p >= input$vol_p, "Significant", "NS")
    df$direction <- ifelse(df$log2FC > 0, "Upregulated", "Downregulated")
    df$color_group <- ifelse(df$sig == "NS", "NS", df$direction)
    
    hl_genes <- find_genes(strsplit(input$vol_hl, "\n")[[1]], a5_genes, a5_genes_upper)
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
    
    if (!is.null(vol_rv$mode) && vol_rv$mode == "paired") {
      color_map <- c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "NS" = "grey70")
      dir_labels <- c("Up in Metastasis", "Up in Primary", "Not Significant")
    } else {
      color_map <- c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "NS" = "grey70")
      dir_labels <- c("Upregulated", "Downregulated", "Not Significant")
    }
    
    p <- ggplot(df, aes(log2FC, neglog10p)) +
      geom_point(aes(color = color_group), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = color_map, name = "Direction",
                         breaks = c("Upregulated", "Downregulated", "NS"),
                         labels = dir_labels) +
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
    filename = function() "volcano_A5.pdf",
    content = function(file) { pdf(file, 10, 8); print(vol_plot_obj()); dev.off() }
  )
  
  output$vol_dl_data <- downloadHandler(
    filename = function() "DEG_results_A5.csv",
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
    df <- a5_meta
    
    if (input$surv_ep == "os") {
      df$time <- as.numeric(df$OS_MONTHS)
      df$status <- as.numeric(df$OS_EVENT)
    } else {
      df$time <- as.numeric(df$MFS_TIME)
      df$status <- as.numeric(df$MFS_EVENT)
    }
    
    df <- df[!is.na(df$time) & df$time > 0 & !is.na(df$status), ]
    
    if (input$surv_by == "gene") {
      req(input$surv_gene %in% a5_genes, input$surv_gene != "")
      df <- df[df$SAMPLE_ID %in% colnames(a5_matrix), ]
      ev <- as.numeric(a5_matrix[input$surv_gene, df$SAMPLE_ID])
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
    filename = function() "survival_A5.pdf",
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
    alterations <- c()
    if (isTRUE(input$onco_tert)) alterations <- c(alterations, "TERT_promoter")
    if (isTRUE(input$onco_atrx)) alterations <- c(alterations, "ATRX_loss")
    if (isTRUE(input$onco_chromo)) alterations <- c(alterations, "Chromothripsis")
    
    if (length(alterations) == 0) {
      plot.new()
      text(0.5, 0.5, "Select at least one alteration to display", cex = 1.5)
      return()
    }
    
    mat <- matrix("", nrow = length(alterations), ncol = nrow(a5_meta))
    rownames(mat) <- alterations; colnames(mat) <- a5_meta$SAMPLE_ID
    
    for (i in seq_len(nrow(a5_meta))) {
      if ("TERT_promoter" %in% alterations && !is.na(a5_meta$TERT_ATRX_STATUS[i]) && a5_meta$TERT_ATRX_STATUS[i] == "TERT_altered")
        mat["TERT_promoter", a5_meta$SAMPLE_ID[i]] <- "Altered"
      if ("ATRX_loss" %in% alterations && !is.na(a5_meta$TERT_ATRX_STATUS[i]) && a5_meta$TERT_ATRX_STATUS[i] == "ATRX_altered")
        mat["ATRX_loss", a5_meta$SAMPLE_ID[i]] <- "Altered"
      if ("Chromothripsis" %in% alterations && !is.na(a5_meta$CHROMOTHRIPSIS[i]) && a5_meta$CHROMOTHRIPSIS[i] == "Yes")
        mat["Chromothripsis", a5_meta$SAMPLE_ID[i]] <- "Altered"
    }
    
    mat <- mat[, order(colSums(mat != ""), decreasing = TRUE), drop = FALSE]
    freq <- round(100 * rowSums(mat != "") / ncol(mat), 1)
    
    annot_vars <- c()
    if (isTRUE(input$onco_lesion)) annot_vars <- c(annot_vars, "LESION_TYPE")
    if (isTRUE(input$onco_location)) annot_vars <- c(annot_vars, "LOCATION_GROUP")
    if (isTRUE(input$onco_mutation)) annot_vars <- c(annot_vars, "MUTATION_TYPE")
    if (isTRUE(input$onco_htn)) annot_vars <- c(annot_vars, "HTN")
    if (isTRUE(input$onco_sex)) annot_vars <- c(annot_vars, "SEX")
    if (isTRUE(input$onco_size)) annot_vars <- c(annot_vars, "TUMOR_SIZE")
    if (isTRUE(input$onco_ki67)) annot_vars <- c(annot_vars, "KI67_GROUP")
    if (isTRUE(input$onco_grantham)) annot_vars <- c(annot_vars, "GRANTHAM_CLASS")
    if (isTRUE(input$onco_prior)) annot_vars <- c(annot_vars, "PRIOR_TREATMENT")
    
    ann_cols <- list(LESION_TYPE = a5_lesion_colors, LOCATION_GROUP = a5_location_colors,
                     MUTATION_TYPE = a5_mutation_colors, HTN = a5_binary_colors,
                     SEX = a5_sex_colors, TUMOR_SIZE = a5_tumor_size_colors,
                     KI67_GROUP = a5_ki67_colors, GRANTHAM_CLASS = a5_grantham_colors,
                     PRIOR_TREATMENT = a5_prior_treatment_colors)
    
    ha <- NULL
    if (length(annot_vars) > 0) {
      ann_df <- a5_meta[match(colnames(mat), a5_meta$SAMPLE_ID), annot_vars, drop = FALSE]
      rownames(ann_df) <- colnames(mat)
      ann_cols_filtered <- ann_cols[intersect(names(ann_cols), annot_vars)]
      ha <- HeatmapAnnotation(df = ann_df, col = ann_cols_filtered, na_col = "grey90",
                              annotation_name_side = "left")
    }
    
    alter_fun <- list(background = alter_graphic("rect", fill = "#EEEEEE"),
                      Altered = alter_graphic("rect", fill = "#E41A1C"))
    
    oncoPrint(mat, alter_fun = alter_fun, col = c(Altered = "#E41A1C"),
              show_column_names = FALSE, alter_fun_is_vectorized = FALSE,
              top_annotation = ha,
              right_annotation = rowAnnotation(Freq = anno_barplot(freq, width = unit(2.5, "cm"))),
              column_title = paste0("A5 SDHB Genomic Alterations (", ncol(mat), " samples)"))
  })
  
  output$onco_plot <- renderPlot({ onco_plot_obj() })
  
  output$onco_table <- renderDT({
    ct <- data.frame(
      Alteration = c("TERT promoter", "ATRX loss", "Chromothripsis"),
      Count = c(sum(a5_meta$TERT_ATRX_STATUS == "TERT_altered", na.rm = TRUE),
                sum(a5_meta$TERT_ATRX_STATUS == "ATRX_altered", na.rm = TRUE),
                sum(a5_meta$CHROMOTHRIPSIS == "Yes", na.rm = TRUE))
    )
    ct$Pct <- paste0(round(100 * ct$Count / nrow(a5_meta), 1), "%")
    datatable(ct, options = list(pageLength = 10, dom = "t"), rownames = FALSE)
  })
  
  output$onco_dl <- downloadHandler(
    filename = function() "oncoprint_A5.pdf",
    content = function(file) { pdf(file, 14, 8); draw(onco_plot_obj()); dev.off() }
  )

  # ==========================================================================
  # PAIRED ANALYSIS
  # ==========================================================================
  
  output$paired_info_table <- renderDT({
    info_df <- data.frame(
      Patient = names(a5_paired_info),
      Primary = sapply(a5_paired_info, function(x) paste(x$primary, collapse = ", ")),
      Metastasis = sapply(a5_paired_info, function(x) paste(x$metastasis, collapse = ", ")),
      SDHB_Variant = sapply(a5_paired_info, function(x) x$variant)
    )
    datatable(info_df, options = list(pageLength = 5, dom = "t"), rownames = FALSE)
  })
  
  paired_plot_obj <- reactive({
    req(input$paired_gene %in% a5_genes, input$paired_gene != "")
    
    plot_data <- data.frame()
    line_data <- data.frame()
    patients_to_show <- if (input$paired_show_all) names(a5_paired_info) else input$paired_patient
    
    for (pt in patients_to_show) {
      info <- a5_paired_info[[pt]]
      for (p_sample in info$primary) {
        if (p_sample %in% colnames(a5_matrix)) {
          plot_data <- rbind(plot_data, data.frame(
            Patient = pt, Sample = p_sample, Type = "Primary",
            Expression = a5_matrix[input$paired_gene, p_sample], Variant = info$variant))
        }
      }
      for (m_sample in info$metastasis) {
        if (m_sample %in% colnames(a5_matrix)) {
          plot_data <- rbind(plot_data, data.frame(
            Patient = pt, Sample = m_sample, Type = "Metastasis",
            Expression = a5_matrix[input$paired_gene, m_sample], Variant = info$variant))
        }
      }
      for (p_sample in info$primary) {
        for (m_sample in info$metastasis) {
          if (p_sample %in% colnames(a5_matrix) && m_sample %in% colnames(a5_matrix)) {
            line_data <- rbind(line_data, data.frame(
              Patient = pt, x_start = "Primary", x_end = "Metastasis",
              y_start = a5_matrix[input$paired_gene, p_sample],
              y_end = a5_matrix[input$paired_gene, m_sample]))
          }
        }
      }
    }
    
    req(nrow(plot_data) > 0)
    plot_data$Type <- factor(plot_data$Type, levels = c("Primary", "Metastasis"))
    patient_colors <- c("E143"="#E41A1C", "E146"="#377EB8", "E158"="#4DAF4A", "E159"="#984EA3", "E225"="#FF7F00")
    
    ggplot() +
      geom_segment(data = line_data, aes(x = x_start, xend = x_end, y = y_start, yend = y_end, color = Patient),
                   linewidth = 1, alpha = 0.7) +
      geom_point(data = plot_data, aes(x = Type, y = Expression, color = Patient, shape = Patient), size = 4, alpha = 0.9) +
      scale_color_manual(values = patient_colors) +
      theme_classic(base_size = 14) +
      labs(x = NULL, y = paste0(input$paired_gene, " Expression (log2 CPM)"),
           title = paste0(input$paired_gene, " - Primary vs Metastasis"),
           subtitle = "Lines connect samples from the same patient") +
      theme(legend.position = "right", axis.text.x = element_text(size = 14, face = "bold"))
  })
  
  output$paired_plot <- renderPlot({ req(input$paired_go > 0 || input$paired_gene != ""); paired_plot_obj() })
  
  output$paired_stats <- renderPrint({
    req(input$paired_gene %in% a5_genes, input$paired_gene != "")
    primary_vals <- c(); met_vals <- c()
    for (pt in names(a5_paired_info)) {
      info <- a5_paired_info[[pt]]
      for (p_sample in info$primary) if (p_sample %in% colnames(a5_matrix)) primary_vals <- c(primary_vals, a5_matrix[input$paired_gene, p_sample])
      for (m_sample in info$metastasis) if (m_sample %in% colnames(a5_matrix)) met_vals <- c(met_vals, a5_matrix[input$paired_gene, m_sample])
    }
    cat("Primary (n=", length(primary_vals), "): mean=", round(mean(primary_vals), 2), ", median=", round(median(primary_vals), 2), "\n", sep = "")
    cat("Metastasis (n=", length(met_vals), "): mean=", round(mean(met_vals), 2), ", median=", round(median(met_vals), 2), "\n", sep = "")
    if (length(primary_vals) >= 3 && length(met_vals) >= 3) {
      wt <- wilcox.test(primary_vals, met_vals)
      cat("\nWilcoxon p =", format(wt$p.value, digits = 3), "\n")
    }
  })
  
  output$paired_dl <- downloadHandler(
    filename = function() paste0("paired_", input$paired_gene, ".pdf"),
    content = function(file) { pdf(file, 10, 7); print(paired_plot_obj()); dev.off() }
  )
  
  output$paired_dl_csv <- downloadHandler(
    filename = function() paste0("paired_analysis_", input$paired_gene, ".csv"),
    content = function(file) {
      plot_data <- data.frame()
      for (pt in names(a5_paired_info)) {
        info <- a5_paired_info[[pt]]
        for (p_sample in info$primary) {
          if (p_sample %in% colnames(a5_matrix)) {
            plot_data <- rbind(plot_data, data.frame(
              Patient = pt, Sample = p_sample, Type = "Primary",
              Expression = a5_matrix[input$paired_gene, p_sample], SDHB_Variant = info$variant))
          }
        }
        for (m_sample in info$metastasis) {
          if (m_sample %in% colnames(a5_matrix)) {
            plot_data <- rbind(plot_data, data.frame(
              Patient = pt, Sample = m_sample, Type = "Metastasis",
              Expression = a5_matrix[input$paired_gene, m_sample], SDHB_Variant = info$variant))
          }
        }
      }
      write.csv(plot_data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
