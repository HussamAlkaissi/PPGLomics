# PPGLomics Home
# Landing page for PPGLomics Explorer
# Developed by Hussam Alkaissi, MD MS
# With assistance from Claude AI

library(shiny)

ui <- navbarPage(
  title = tags$span(
    tags$strong("PPGLomics", style = "color: white;"),
    tags$span(" | ", style = "color: #888; margin: 0 5px;"),
    tags$span("Alkaissi Lab", style = "color: #E8654B; font-weight: 500;")
  ),
  id = "main_nav",
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  tags$head(tags$style(HTML("
    .main-container {
      max-width: 1000px;
      margin: 0 auto;
      padding: 40px 20px;
    }
    .title-section {
      text-align: center;
      margin-bottom: 50px;
    }
    .main-title {
      color: #2C3E50;
      font-size: 3em;
      font-weight: bold;
      margin-bottom: 10px;
    }
    .subtitle {
      color: #666;
      font-size: 1.2em;
      margin-bottom: 20px;
    }
    .cards-container {
      display: flex;
      justify-content: center;
      gap: 40px;
      flex-wrap: wrap;
      margin-bottom: 50px;
    }
    .dataset-card {
      background: #ffffff;
      border-radius: 15px;
      padding: 35px;
      width: 380px;
      text-align: center;
      box-shadow: 0 5px 20px rgba(0,0,0,0.1);
      transition: transform 0.3s ease, box-shadow 0.3s ease;
      border: 1px solid #eee;
    }
    .dataset-card:hover {
      transform: translateY(-5px);
      box-shadow: 0 10px 30px rgba(0,0,0,0.15);
    }
    .card-tcga {
      border-top: 5px solid #3498db;
    }
    .card-a5 {
      border-top: 5px solid #e74c3c;
    }
    .card-title {
      font-size: 1.6em;
      font-weight: bold;
      margin-bottom: 10px;
    }
    .card-tcga .card-title { color: #3498db; }
    .card-a5 .card-title { color: #e74c3c; }
    .card-samples {
      font-size: 2.2em;
      font-weight: bold;
      color: #2c3e50;
      margin-bottom: 5px;
    }
    .card-desc {
      color: #666;
      font-size: 0.95em;
      margin-bottom: 15px;
      line-height: 1.5;
    }
    .card-features {
      text-align: left;
      margin-bottom: 20px;
      padding: 12px 15px;
      background: #f8f9fa;
      border-radius: 8px;
    }
    .card-features li {
      color: #555;
      margin-bottom: 6px;
      font-size: 0.9em;
    }
    .btn-explore {
      display: inline-block;
      padding: 12px 35px;
      font-size: 1em;
      font-weight: bold;
      text-decoration: none;
      border-radius: 25px;
      transition: all 0.3s ease;
      border: none;
      cursor: pointer;
    }
    .btn-tcga {
      background: #3498db;
      color: white !important;
    }
    .btn-tcga:hover {
      background: #2980b9;
      text-decoration: none;
    }
    .btn-a5 {
      background: #e74c3c;
      color: white !important;
    }
    .btn-a5:hover {
      background: #c0392b;
      text-decoration: none;
    }
    .reference {
      font-size: 0.85em;
      color: #888;
      margin-top: 10px;
      margin-bottom: 15px;
    }
    .footer {
      text-align: center;
      color: #666;
      margin-top: 30px;
      padding-top: 20px;
      border-top: 1px solid #eee;
    }
  "))),
  
  tabPanel("",
    div(class = "main-container",
      div(class = "title-section",
        h1(class = "main-title", "PPGLomics Explorer"),
        p(class = "subtitle", "Integrated Transcriptomics Analysis for Pheochromocytoma and Paraganglioma")
      ),
      
      div(class = "cards-container",
        # TCGA Card
        div(class = "dataset-card card-tcga",
          h2(class = "card-title", "TCGA-PCPG"),
          div(class = "card-samples", "160"),
          p(style = "color: #666; margin-bottom: 12px;", "samples"),
          p(class = "card-desc", "Multiple genotypes including SDHx, VHL, RET, NF1, HRAS, and more. mRNA clusters and methylation patterns."),
          div(class = "card-features",
            tags$ul(
              tags$li("18 driver genotypes"),
              tags$li("4 mRNA expression clusters"),
              tags$li("Methylation subtypes"),
              tags$li("Survival analysis")
            )
          ),
          p(class = "reference", "Fishbein et al., Cancer Cell 2017"),
          tags$a(href = "https://alkaissilab.shinyapps.io/PPGLomics_TCGA/", 
                 class = "btn-explore btn-tcga", target = "_blank", "Explore TCGA-PCPG →")
        ),
        
        # A5 Card
        div(class = "dataset-card card-a5",
          h2(class = "card-title", "A5 Consortium"),
          div(class = "card-samples", "91"),
          p(style = "color: #666; margin-bottom: 12px;", "samples"),
          p(class = "card-desc", "SDHB germline mutations only. Includes primary/metastatic pairs, TERT/ATRX status, and Grantham scores."),
          div(class = "card-features",
            tags$ul(
              tags$li("SDHB-deficient tumors"),
              tags$li("Paired primary/metastasis"),
              tags$li("TERT/ATRX alterations"),
              tags$li("Prior treatment status")
            )
          ),
          p(class = "reference", "Flynn et al., Nat Commun 2025"),
          tags$a(href = "https://alkaissilab.shinyapps.io/PPGLomics_A5/", 
                 class = "btn-explore btn-a5", target = "_blank", "Explore A5 →")
        )
      ),
      
      div(class = "footer",
        p(tags$strong("Developed by Hussam Alkaissi, MD MS")),
        p("Clinician-scientist | Pseudohypoxia & HIF signaling research"),
        p(em("With assistance from Claude AI")),
        br(),
        p(style = "font-size: 0.9em; color: #999;", "PPGLomics v1.0 | Alkaissi Lab")
      )
    )
  )
)

server <- function(input, output, session) {
  # No server logic needed for static landing page
}

shinyApp(ui, server)
