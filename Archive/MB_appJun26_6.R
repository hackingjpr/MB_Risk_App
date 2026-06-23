source("./MB_app_source_functionsJun26.R")

library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(fs)
library(later)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  useAttendant(),
  
  dashboardPage(
    skin = "blue",
    
    dashboardHeader(
      title = span("MB Risk Calculator", style = "font-weight:600;")
    ),
    
    dashboardSidebar(
      div(style = "text-align:center; padding:10px;",
          img(src='Dial_brain.png', width = "80%"),
          tags$hr()
      ),
      
      h4("Step 1: Upload Data"),
      fileInput("idatFile", NULL, multiple = TRUE, accept = ".idat"),
      
      h4("Step 2: Select Model"),
      selectInput("metagenes", NULL,
                  c("SHH", "G3_G4_no_sub", "G3_G4_sub")),
      
      h4("Step 3: Run Analysis"),
      
      actionBttn("bttn1",
                 label = tagList(icon("play"), "Run Analysis"),
                 color = "primary",
                 style = "fill"),
      
      br(), br(),
      attendantBar("progress-bar"),
      br(),
      
      actionBttn("bttn2",
                 label = tagList(icon("redo"), "Reset"),
                 color = "danger",
                 style = "bordered",
                 size = "sm")
    ),
    
    dashboardBody(
      tabsetPanel(id = "tabs",
                  
                  tabPanel("Info",
                           fluidRow(
                             box(width = 12, title = "Introduction",
                                 includeMarkdown("./introduction.md"))
                           )
                  ),
                  
                  tabPanel("Tutorial",
                           fluidRow(
                             box(width = 12, title = "Tutorial",
                                 includeMarkdown("./Tutorial/tutorial.md"))
                           )
                  ),
                  
                  tabPanel("Results",
                           
                           fluidRow(
                             box(width = 4, title = "Metagene",
                                 h3(textOutput("metagenechoice"))),
                             
                             box(width = 4, title = "Sample",
                                 h3(textOutput("sample"))),
                             
                             box(width = 4, title = "Interpretation",
                                 textOutput("interpretation"))
                           ),
                           
                           fluidRow(
                             box(width = 12, title = "Clinical Inputs",
                                 uiOutput("sample_clinical_inputs"))
                           ),
                           
                           fluidRow(
                             box(width = 4, title = "Risk Distribution",
                                 plotOutput("figure")),
                             
                             box(width = 4, title = "Survival Curve",
                                 plotOutput("survivalPlot")),
                             
                             box(width = 4, title = "Risk Values",
                                 DTOutput("Mval"))
                           )
                  )
      )
    )
  )
)

# =========================
# SERVER
# =========================

server <- function(input, output, session) {
  
  hideTab("tabs", "Results")
  
  att <- Attendant$new("progress-bar")
  
  dir.create("./temp", showWarnings = FALSE)
  
  Mvals <- reactiveValues(data = NULL)
  
  
  observeEvent(input$bttn1, {
    
    req(input$idatFile)
    
    # ✅ CAPTURE INPUTS FIRST (VERY IMPORTANT)
    idat_files <- input$idatFile
    metagene_choice <- input$metagenes
    
    tempDIR <- fs::path("./temp",
                        paste0(sample(letters,10,TRUE), collapse=""))
    dir.create(tempDIR)
    
    att$set(5, "Starting...")
    
    later(function() {
      
      # ✅ USE CAPTURED VALUES (NOT input$...)
      file.copy(idat_files$datapath,
                file.path(tempDIR, idat_files$name))
      
      att$set(20, "Files uploaded")
      
      temp.base <- get_basenames(tempDIR)
      
      att$set(40, "Processing data")
      
      temp.processed <- process_idats(temp.base)
      
      meta <- switch(metagene_choice,
                     "SHH" = SHH,
                     "G3_G4_sub" = G3_G4_sub,
                     "G3_G4_no_sub" = G3_G4_no_sub)
      
      test.res <- extract.metagene(
        as.character(meta[[1]]$genes),
        as.numeric(meta[[1]]$weights),
        beta2m(temp.processed$betas),
        as.numeric(meta[[2]])
      )
      
      test.res <- round(test.res, 3)
      Mvals$data <- test.res   # ✅ reactiveValues is fine
      
      # ✅ Render outputs OUTSIDE reactive context is OK
      output$Mval <- renderDT({
        datatable(
          Mvals$data,
          selection = list(mode = "single", target = "row")
        )
      })
      
      att$set(100, "Complete")
      
    }, delay = 0.1)
  })
  
  observeEvent(input$bttn2, {
    session$reload()
  })
}

shinyApp(ui, server)