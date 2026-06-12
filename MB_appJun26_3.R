source("./MB_app_source_functionsJun26.R")

library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(gridExtra)
library(grid)
library(fs)

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
  
  # ===== MAIN TABS =====
  navset_tab(
    id = "tabs",
    
    # ✅ DEFAULT LANDING PAGE
    nav("Introduction",
        card(
          card_header("Overview"),
          includeMarkdown("./introduction.md")
        )
    ),
    
    # ✅ ALWAYS AVAILABLE
    nav("Tutorial",
        card(
          card_header("How to Use"),
          includeMarkdown("./Tutorial/tutorial.md")
        )
    ),
    
    # ✅ RESULTS (hidden initially)
    nav("Results",
        value = "results",
        
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
    
    # ✅ DOWNLOAD TAB
    nav("Download",
        value = "download",
        
        card(
          card_header("Export Results"),
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
  
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  
  att <- Attendant$new("progress-bar")
  
  dir.create("./temp", showWarnings = FALSE)
  tempDIR <- fs::path("./temp", paste0(sample(letters, 10, TRUE), collapse = ""))
  dir.create(tempDIR)
  
  Mvals <- reactiveValues(data = NULL)
  report_text <- reactiveVal("")
  
  observeEvent(input$bttn1, {
    
    req(input$idatFile)
    
    start_time <- Sys.time()
    
    att$set(0, "Starting...")
    
    # STEP 1
    file.copy(input$idatFile$datapath,
              file.path(tempDIR, input$idatFile$name))
    
    att$set(15, "Files uploaded")
    
    # STEP 2
    temp.base <- get_basenames(tempDIR)
    att$set(25, "Preparing data")
    
    # STEP 3 (heavy)
    temp.processed <- process_idats(temp.base)
    
    elapsed <- round(as.numeric(Sys.time() - start_time, units = "secs"))
    att$set(60, paste0("Processing methylation data (", elapsed, "s elapsed)"))
    
    # STEP 4
    meta <- switch(input$metagenes,
                   "SHH" = SHH,
                   "G3_G4_sub" = G3_G4_sub,
                   "G3_G4_no_sub" = G3_G4_no_sub)
    
    att$set(70, "Applying metagene model")
    
    # STEP 5
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
    
    # AUTO REPORT TEXT
    report_text(paste(
      "Medulloblastoma Risk Analysis Report\n\n",
      "Metagene Model:", input$metagenes, "\n\n",
      "Summary:\n",
      interpretation, "\n\n",
      "The analysis was performed using methylation array data processed ",
      "through a predefined metagene model.\n\n",
      "Average risk score:", round(risk, 3), "\n\n",
      "This tool is intended for research use only and should not be used ",
      "for clinical decision making."
    ))
    
    att$set(100, "Complete")
    
    output$metagenechoice <- renderText(input$metagenes)
    output$sample <- renderText(input$Mval_cell_clicked$value)
  })
  
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = TRUE)
    session$reload()
  })
  
  # =========================
  # DOWNLOAD HANDLER
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
        
        # TITLE
        grid.text("Medulloblastoma Risk Report",
                  x = 0.5, y = 0.95,
                  gp = gpar(fontsize = 18, fontface = "bold"))
        
        # REPORT TEXT
        grid.text(report_text(),
                  x = 0.05, y = 0.80,
                  just = "left",
                  gp = gpar(fontsize = 11))
        
        # TABLE
        grid.newpage()
        gridExtra::grid.table(round(Mvals$data, 3))
        
        # PLOT
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