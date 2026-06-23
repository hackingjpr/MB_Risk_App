source("./MB_app_source_functionsJun26.R")

library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(gridExtra)
library(grid)
library(fs)
library(future)
library(promises)

# Set max file size limit
options(shiny.maxRequestSize = 200 * 1024^2)

# Enable async processing
plan(multisession)

# =========================
# UI
# =========================

ui <- page_sidebar(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    base_font = font_google("Inter")
  ),
  
  useAttendant(),
  
  sidebar = sidebar(
    
    div(style = "text-align:center;",
        img(src='MB_risk.jpg', width = "70%")),
    
    hr(),
    
    h5("Step 1: Upload Data"),
    fileInput("idatFile", "Upload IDAT files", multiple = TRUE),
    
    h5("Step 2: Select Model"),
    selectInput("metagenes", "Metagene Set",
                c("SHH", "G3_G4_sub", "G3_G4_no_sub")),
    
    h5("Step 3: Run Analysis"),
    actionButton("bttn1", "Run Analysis", class = "btn-primary"),
    
    br(), br(),
    
    attendantBar("progress-bar"),
    
    hr(),
    
    actionButton("bttn2", "Reset", class = "btn-outline-danger btn-sm")
  ),
  
  navset_tab(
    id = "tabs",
    
    nav("Introduction",
        card(includeMarkdown("./introduction.md"))
    ),
    
    nav("Tutorial",
        card(includeMarkdown("./Tutorial/tutorial.md"))
    ),
    
    nav("Results", value = "results",
        
        card(
          card_header("Risk Table"),
          DTOutput("Mval")
        ),
        
        card(
          card_header("Risk Plot"),
          plotOutput("figure", height = "450px")
        ),
        
        div(class = "row",
            div(class = "col-md-6",
                card(
                  card_header("Selection Info"),
                  textOutput("metagenechoice"),
                  textOutput("sample")
                )
            ),
            div(class = "col-md-6",
                card(
                  card_header("Clinical Interpretation"),
                  textOutput("interpretation")
                )
            )
        )
    ),
    
    nav("Download", value = "download",
        card(
          textInput("filename", "Filename", "MB_report"),
          radioButtons("download", "File Type", c("csv", "pdf")),
          downloadButton("down", "Download")
        )
    )
  )
)

# =========================
# SERVER
# =========================

server <- function(input, output, session) {
  
  att <- Attendant$new("progress-bar")
  
  # Hide result tabs initially
  hideTab("tabs", "results")
  hideTab("tabs", "download")
  
  dir.create("./temp", showWarnings = FALSE)
  
  Mvals <- reactiveValues(data = NULL)
  report_text <- reactiveVal("")
  
  observeEvent(input$bttn1, {
    
    req(input$idatFile)
    
    tempDIR <- fs::path("./temp",
                        paste0(sample(letters, 10, TRUE), collapse = ""))
    dir.create(tempDIR)
    
    att$set(0, "Starting...")
    
    # Copy files (fast, keep synchronous)
    file.copy(input$idatFile$datapath,
              file.path(tempDIR, input$idatFile$name))
    
    att$set(10, "Files uploaded")
    
    temp.base <- get_basenames(tempDIR)
    
    att$set(20, "Preparing analysis...")
    
    # ✅ ASYNC BLOCK STARTS
    future({
      
      # HEAVY STEP runs in background
      temp.processed <- process_idats(temp.base)
      
      list(temp.processed = temp.processed)
      
    }) %...>% (function(res) {
      
      att$set(60, "Processing complete")
      
      # MODEL SELECTION
      meta <- switch(input$metagenes,
                     "SHH" = SHH,
                     "G3_G4_sub" = G3_G4_sub,
                     "G3_G4_no_sub" = G3_G4_no_sub)
      
      att$set(70, "Applying model")
      
      # CALCULATION
      test.res <- extract.metagene(
        as.character(meta[[1]]$genes),
        as.numeric(meta[[1]]$weights),
        beta2m(res$temp.processed$betas),
        as.numeric(meta[[2]])
      )
      
      test.res <- round(test.res, 3)
      Mvals$data <- test.res
      
      # TABLE
      output$Mval <- renderDT({
        datatable(Mvals$data, options = list(pageLength = 10))
      })
      
      # PLOT
      figure.input <- test.res$Risk_Value
      names(figure.input) <- rownames(test.res)
      
      output$figure <- renderPlot({
        if (input$metagenes == "SHH") {
          generate_figure_highlight_SHH(figure.input, 1)
        } else if (input$metagenes == "G3_G4_sub") {
          generate_figure_highlight_G3_G4_sub(figure.input, 1)
        } else {
          generate_figure_highlight_G3_G4_no_sub(figure.input, 1)
        }
      })
      
      # INTERPRETATION
      risk <- mean(test.res$Risk_Value, na.rm = TRUE)
      
      interpretation <- if (risk < 0.3) {
        "Low-risk methylation profile. Associated with favourable prognosis."
      } else if (risk < 0.7) {
        "Intermediate-risk profile. Clinical context required."
      } else {
        "High-risk profile. Associated with poorer outcomes."
      }
      
      output$interpretation <- renderText(interpretation)
      
      report_text(paste(
        "Medulloblastoma Risk Analysis Report\n\n",
        "Metagene Model:", input$metagenes, "\n\n",
        "Summary:\n", interpretation, "\n\n",
        "Mean risk score:", round(risk, 3)
      ))
      
      att$set(100, "Complete")
      
      output$metagenechoice <- renderText(input$metagenes)
      output$sample <- renderText(input$Mval_cell_clicked$value)
      
      # ✅ AUTO SWITCH TO RESULTS
      showTab("tabs", "results")
      showTab("tabs", "download")
      
      nav_select(id = "tabs", selected = "results", session = session)
      
    }) %...!% (function(e) {
      
      # ERROR HANDLING
      att$set(100, "Error occurred")
      showNotification(paste("Error:", e$message), type = "error")
      
    })
  })
  
  observeEvent(input$bttn2, {
    session$reload()
  })
  
  # =========================
  # DOWNLOAD
  # =========================
  
  output$down <- downloadHandler(
    
    filename = function() {
      paste(input$filename, Sys.Date(),
            ifelse(input$download == "csv", ".csv", ".pdf"),
            sep = "_")
    },
    
    content = function(file) {
      
      if (input$download == "csv") {
        
        write.csv(Mvals$data, file, row.names = TRUE)
        
      } else {
        
        pdf(file, width = 11, height = 8.5)
        
        grid.text("Medulloblastoma Risk Report",
                  x = 0.5, y = 0.95,
                  gp = gpar(fontsize = 18, fontface = "bold"))
        
        grid.text(report_text(),
                  x = 0.05, y = 0.85,
                  just = "left",
                  gp = gpar(fontsize = 11))
        
        grid.newpage()
        grid.table(round(Mvals$data, 3))
        
        grid.newpage()
        
        figure.input <- Mvals$data$Risk_Value
        names(figure.input) <- rownames(Mvals$data)
        
        if (input$metagenes == "SHH") {
          print(generate_figure_highlight_SHH(figure.input, 1))
        } else if (input$metagenes == "G3_G4_sub") {
          print(generate_figure_highlight_G3_G4_sub(figure.input, 1))
        } else {
          print(generate_figure_highlight_G3_G4_no_sub(figure.input, 1))
        }
        
        dev.off()
      }
    }
  )
}

shinyApp(ui, server)
