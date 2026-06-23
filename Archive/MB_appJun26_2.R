source("./MB_app_source_functionsJun26.R")

library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(fs)

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
                 size = "sm"),
      
      tags$hr(),
      
      div(style = "font-size: 12px;",
          HTML("
            <b>Author:</b> James Hacking<br/>
            <b>Date:</b> 12-06-2026<br/>
            <b>Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic or clinical use.
          "),
          HTML(
            paste(
              h5(
                a("Email: james.hacking@ncl.ac.uk", href = "mailto:james.hacking@ncl.ac.uk")
              )
              
            ))
      )
    ),
    
    dashboardBody(
      tags$head(
        
        # Load Inter font
        tags$link(
          href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600&display=swap",
          rel = "stylesheet"
        ),
        
        # Custom styles
        tags$style(HTML("
      body {
        font-family: 'Inter', 'Segoe UI', Arial;
        font-weight: 400;
      }

      .content-wrapper {
        background-color: #f5f7fa;
      }

      .box {
        border-radius: 10px;
      }

      /* Markdown styling */
      .box-body {
        font-size: 14px;
        line-height: 1.5;
        font-weight: 400;
      }

      .box-body h1 { font-size: 20px; font-weight: 600; }
      .box-body h2 { font-size: 18px; font-weight: 600; }
      .box-body h3 { font-size: 16px; font-weight: 500; }

      .box-body p {
        margin-bottom: 8px;
      }
    "))
      ),
      
      tabsetPanel(id = "tabs",
                  
                  tabPanel("Info",
                           fluidRow(
                             box(
                               width = 12,
                               solidHeader = TRUE,
                               status = "primary",
                               title = tagList(icon("info-circle"), "Introduction"),
                               
                               div(style = "padding: 15px; font-size: 15px; line-height: 1.6;",
                                   includeMarkdown("./introduction.md")
                               )
                             )
                           )
                  ),
                  
                  tabPanel("Tutorial",
                           fluidRow(
                             box(
                               width = 12,
                               solidHeader = TRUE,
                               status = "primary",
                               title = tagList(icon("book"), "Tutorial"),
                               
                               div(
                                 style = "
          padding: 15px;
          font-size: 14px;
          line-height: 1.5;
          font-weight: 400;
        ",
                                 includeMarkdown("./Tutorial/tutorial.md")
                               )
                             )
                           )
                  ),
                  
                  tabPanel("Results",
                           
                                    
                                    # ✅ NEW: Results explanation box
                                    fluidRow(
                                      box(
                                        width = 12,
                                        solidHeader = TRUE,
                                        collapsible = TRUE,
                                        collapsed = TRUE,
                                        status = "info",
                                        title = tagList(icon("info-circle"), "Understanding your results"),
                                        
                                        div(
                                          style = "
          padding: 15px;
          font-size: 14px;
          line-height: 1.5;
          font-weight: 400;
        ",
                                          includeMarkdown("./results_blurb.md")
                                        )
                                      )
                                    ),
                                    
                                    # Existing content
                                    fluidRow(
                                      box(
                                        width = 12,
                                        title = tagList(icon("table"), "Risk Values"),
                                        status = "primary", solidHeader = TRUE,
                                        DTOutput("Mval")
                                      )
                                    ),
                           
                           fluidRow(
                             box(
                               title = tagList(icon("chart-line"), "Risk Plot"),
                               status = "warning", solidHeader = TRUE,
                               plotOutput("figure", height = "320px")
                             ),
                             
                             box(
                               title = tagList(icon("info-circle"), "Selections"),
                               status = "success", solidHeader = TRUE,
                               h4("Metagene:"), textOutput("metagenechoice"),
                               h4("Sample:"), textOutput("sample")
                             )
                           )
                  ),
                  
                  tabPanel("Download",
                           textInput("filename", "Filename", "M-values"),
                           radioButtons("download", "File type", c("csv", "pdf")),
                           downloadButton("down", "Download results")
                  )
      )
    )
  )
)

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  
  hideTab("tabs", "Results")
  hideTab("tabs", "Download")
  
  att <- Attendant$new("progress-bar")
  
  unlink("./temp/", recursive = TRUE)
  dir.create("./temp", showWarnings = FALSE)
  
  tempDIR <- fs::path("./temp/", paste0(sample(letters,10,TRUE), collapse=""))
  dir.create(tempDIR)
  
  Mvals <- reactiveValues(data = NULL)
  
  observeEvent(input$bttn1, {
    
    req(input$idatFile)
    
    att$set(0, "Starting...")
    
    # STEP 1: File copy
    in.files <- input$idatFile$datapath
    out.files <- file.path(tempDIR, input$idatFile$name)
    file.copy(in.files, out.files)
    
    att$set(15, "Files uploaded")
    
    # STEP 2: Prepare filenames
    temp.base <- get_basenames(tempDIR)
    att$set(25, "Preparing data")
    
    # STEP 3: Heavy processing
    ptm <- proc.time()
    temp.processed <- process_idats(temp.base)
    
    att$set(60, "Processing methylation data")
    
    # STEP 4: Select model
    meta <- switch(input$metagenes,
                   "SHH" = SHH,
                   "G3_G4_sub" = G3_G4_sub,
                   "G3_G4_no_sub" = G3_G4_no_sub)
    
    att$set(70, "Applying model")
    
    # STEP 5: Calculate risk
    test.res <- extract.metagene(
      as.character(meta[[1]]$genes),
      as.numeric(meta[[1]]$weights),
      beta2m(temp.processed$betas),
      as.numeric(meta[[2]])
    )
    
    test.res <- round(test.res, 3)
    Mvals$data <- test.res
    
    att$set(85, "Generating outputs")
    
    # TABLE
    output$Mval <- renderDT({
      datatable(test.res, options = list(pageLength = 10))
    })
    
    # FIGURE
    figure.input <- test.res$Risk_Value
    names(figure.input) <- rownames(test.res)
    
    att$set(95, "Rendering plot")
    
    output$figure <- renderPlot({
      if (input$metagenes == "SHH") {
        generate_figure_highlight_SHH(figure.input, input$Mval_row_last_clicked)
      } else if (input$metagenes == "G3_G4_sub") {
        generate_figure_highlight_G3_G4_sub(figure.input, input$Mval_row_last_clicked)
      } else {
        generate_figure_highlight_G3_G4_no_sub(figure.input, input$Mval_row_last_clicked)
      }
    })
    
    att$set(100, "Complete")
    
    output$metagenechoice <- renderText(input$metagenes)
    output$sample <- renderText(input$Mval_cell_clicked$value)
    
    showTab("tabs", "Results", select = TRUE)
    showTab("tabs", "Download")
    
    output$time <- renderText(proc.time() - ptm)
  })
  
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = TRUE)
    session$reload()
  })
  
  output$down <- downloadHandler(
    filename = function() {
      paste(input$filename, Sys.time(),
            ifelse(input$download == "csv", ".csv", ".pdf"),
            sep = "_")
    },
    
    content = function(file) {
      
      if (input$download == "csv") {
        write.csv(Mvals$data, file, row.names = TRUE)
        
      } else {
        
        pdf(file, width = 14, height = 8)
        gridExtra::grid.table(round(Mvals$data, 3))
        dev.off()
      }
    }
  )
}

shinyApp(ui, server)