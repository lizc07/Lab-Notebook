library(shiny)
load("data.Rdata")
shinyUI(fluidPage(
  titlePanel(title = NULL,windowTitle = "Discovering TCGA-BRCA"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Analysis of genes and miRNAs in TCGA BRCA", style = "color:purple"),
      h5("STEP0. Select genes and miRNAs", style = "font-weight:bold; color:blue"),
      selectInput("sel_gene", label = "select genes", choices = unique(genes), selected = "HIF1A", multiple = F),
      selectInput("sel_mirna", label = "select miRNAs", choices = unique(mirna), selected = "hsa-miR-153-3p|MIMAT0000439", multiple = F),
      
      radioButtons("analysis", label = h5("STEP1. Select Analysis Type", style = "font-weight:bold; color:blue"), choices = list("Differential expression", "Correlation"), selected = "Differential expression",inline = F),
      radioButtons("category", label = h5("STEP2. Select Groups by category", style = "font-weight:bold; color:blue"), choices = list(Tumor_or_Normal = "N_C_TYPES","TNBC_SUBTYPES","PAM50_SUBTYPES","PATHOLOGICAL_STAGES","METASTASIS_STATUS","LYMPH_NODE_STATUS","MENOPAUSE_STATUS","ER_STATUS","PR_STATUS","HER2_STATUS"), selected = "N_C_TYPES",inline = F),
      sliderInput("heat.height", label = ("Plot height"),
                  min = 500, max = 2000, value = 800,ticks = F),
      submitButton("Submit")
      
    ),
   
    mainPanel(
      tabsetPanel(
        tabPanel("Plots",uiOutput("plot.ui")),
        tabPanel("P-value",tableOutput("testTab"))

      )
    )
  )
))