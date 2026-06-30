source("./source_functions.R")

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
                  c("SHH", "Group3/4 (Early)", "Group3/4 (Late)" )),
      
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
                               hr(),
                               uiOutput("covariate_editor"),
                               br(),
                               uiOutput("update_survival_btn")
                             )
                             
                           ),
                           fluidRow(
                             box(
                               width = 12,
                               title = tagList(icon("heartbeat"), "Survival Estimate"),
                               status = "danger", solidHeader = TRUE,
                               div(
                                 style = "padding: 5px 10px; font-size: 14px;",
                                 strong("Selected sample: "),
                                 textOutput("selected_sample_surv", inline = TRUE),
                                 br(),
                                 strong("Estimated survival: "),
                                 textOutput("selected_surv_value", inline = TRUE)
                               ),
                               hr(),
                               plotOutput("survivalFigure", height = "320px")
                             )
                           )
                  ),
                  
                  tabPanel("Download",
                           textInput("filename", "Filename", "Risk-Values"),
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
  
  tempDIR <- fs::path("./temp/", paste0(sample(letters, 10, TRUE), collapse = ""))
  dir.create(tempDIR)
  
  
  
  survInputs <- reactiveValues(
    score      = NULL,
    covariates = NULL,   # per-sample covariate data frame
    indexRow   = 1
  )
  
  # Render the per-row covariate editor when a row is selected
  output$covariate_editor <- renderUI({
    req(survInputs$covariates)
    
    metagene <- input$metagenes
    
    if (metagene == "Group3/4 (Late)") return(NULL)
    
    tagList(
      h4("Sample Covariates"),
      p(style = "font-size:12px; color:#666;",
        "Edit values directly in the table, then click Update."),
      rHandsontableOutput("covariate_table"),
      br()
    )
  })
  output$covariate_table <- renderRHandsontable({
    req(survInputs$covariates)
    
    covs     <- survInputs$covariates
    metagene <- input$metagenes
    
    # Guard — this model has no covariates
    if (metagene == "Group3/4 (Late)") return(NULL)
    
    if (metagene == "SHH") {
      display_df <- data.frame(
        Sample         = covs$sample,
        MYCN_Amplified = as.logical(covs$mycn),
        stringsAsFactors = FALSE
      )
    } else if (metagene == "Group3/4 (Early)") {
      display_df <- data.frame(
        Sample        = covs$sample,
        MYC_Amplified = as.logical(covs$myc),
        Metastatic    = as.logical(covs$mets),
        stringsAsFactors = FALSE
      )
    } else {
      return(NULL)  # fallback safety
    }
    
    rhandsontable(display_df, rowHeaders = FALSE, stretchH = "all") %>%
      hot_col("Sample", readOnly = TRUE) %>%
      hot_cols(halign = "htCenter") %>%
      hot_table(highlightCol = FALSE, highlightRow = TRUE)
  })
  
  # Save edits back to the per-sample covariate df and replot
  observeEvent(input$bttn_update_surv, {
    req(survInputs$covariates, survInputs$score)
    
    metagene <- input$metagenes
    
    # Read current state of the edited table
    if (!is.null(input$covariate_table)) {
      edited <- hot_to_r(input$covariate_table)
      
      message("=== UPDATE BUTTON CLICKED ===")
      message("Metagene: ", metagene)
      message("Edited table:")
      print(edited)
      
      covs <- survInputs$covariates
      
      message("Covariates BEFORE update:")
      print(covs)
      
      if (metagene == "SHH") {
        covs$mycn <- as.integer(edited$MYCN_Amplified)
        
      } else if (metagene == "Group3/4 (Early)") {
        covs$myc  <- as.integer(edited$MYC_Amplified)
        covs$mets <- as.integer(edited$Metastatic)
      }
      
      message("Covariates AFTER update:")
      print(covs)
      
      survInputs$covariates <- covs
    } else {
      message("covariate_table is NULL — table not being read")
    }
  })
  
  # Replot whenever survInputs changes
  observe({
    req(survInputs$score, survInputs$covariates)
    
    output$selected_sample_surv <- renderText({
      req(survInputs$covariates, survInputs$indexRow)
      covs     <- survInputs$covariates
      indexRow <- survInputs$indexRow
      if (indexRow <= nrow(covs)) covs$sample[indexRow] else "None"
    })
    
    output$selected_surv_value <- renderText({
      req(survInputs$score, survInputs$covariates)
      score    <- survInputs$score
      covs     <- survInputs$covariates
      indexRow <- survInputs$indexRow
      metagene <- input$metagenes
      tryCatch({
        if (metagene == "SHH") {
          fit <- survfit(.shh_fit,
                         newdata = data.frame(
                           fixedMG       = as.numeric(score[indexRow]),
                           ConsensusMYCN = as.numeric(covs$mycn[indexRow])
                         ))
          val <- as.numeric(summary(fit, time = 5)$surv)
          paste0(round(val * 100, 1), "% at 5 years")
          
        } else if (metagene == "Group3/4 (Early)") {
          new_row <- data.frame(
            fixedMG      = as.numeric(score[indexRow]),
            M._versus_M. = as.numeric(covs$mets[indexRow]),
            ConsensusMYC = factor(
              paste0("ConsensusMYC=", covs$myc[indexRow]),
              levels = .g34early_fit$xlevels[["strata(ConsensusMYC)"]]
            )
          )
          fit <- survfit(.g34early_fit, newdata = new_row)
          val <- as.numeric(summary(fit, time = 5)$surv)
          paste0(round(val * 100, 1), "% at 5 years")
          
        } else {
          fit <- survfit(.g34late_fit,
                         newdata = data.frame(
                           fixedMG = as.numeric(score[indexRow])
                         ))
          val <- as.numeric(summary(fit, time = 10)$surv)
          paste0(round(val * 100, 1), "% at 10 years")
        }
      }, error = function(e) "Unable to compute")
    })
    
    output$survivalFigure <- renderPlot({
      # Read directly from survInputs here — not captured local variables
      req(survInputs$score, survInputs$covariates)
      score    <- survInputs$score
      covs     <- survInputs$covariates
      indexRow <- survInputs$indexRow
      metagene <- input$metagenes
      
      render_survival_plot(
        metagene = metagene,
        score    = score,
        covs     = covs,
        indexRow = indexRow
      )
    })
  })
  
  Mvals <- reactiveValues(data = NULL)
  
  # ── Observers that must live at the top level ──────────────────────────────
  
  observeEvent(input$Mval_row_last_clicked, {
    survInputs$indexRow <- input$Mval_row_last_clicked
  })
  
  # observeEvent(input$bttn_update_surv, {
  #   req(survInputs$score)
  #   
  # })
  
  
  # ── Main analysis button ───────────────────────────────────────────────────
  
  observeEvent(input$bttn1, {
    
    req(input$idatFile)
    att$set(0, "Starting...")
    
    in.files  <- input$idatFile$datapath
    out.files <- file.path(tempDIR, input$idatFile$name)
    file.copy(in.files, out.files)
    att$set(15, "Files uploaded")
    
    temp.base <- get_basenames(tempDIR)
    att$set(25, "Preparing data")
    
    ptm <- proc.time()
    temp.processed <- process_idats(temp.base)
    att$set(60, "Processing methylation data")
    
    meta <- switch(input$metagenes,
                   "SHH"          = SHH,
                   "Group3/4 (Early)"    = G3_G4_sub,
                   "Group3/4 (Late)" = G3_G4_no_sub)
    att$set(70, "Applying model")
    
    message("=== DEBUG ===")
    message("Metagene selected: ", input$metagenes)
    message("Number of CpGs in meta: ", nrow(meta[[1]]))
    message("Scaling: ", meta[[2]])
    
    test.res <- extract.metagene(
      as.character(meta[[1]]$genes),
      as.numeric(meta[[1]]$weights),
      beta2m(temp.processed$betas),
      as.numeric(meta[[2]])
    )
    test.res   <- round(test.res, 3)
    Mvals$data <- test.res
    
    
    Mvals$metagene <- input$metagenes          # snapshot at compute time 
    survInputs$metagene <- input$metagenes
    
    # Define figure.input BEFORE using it
    figure.input        <- test.res$Risk_Value
    names(figure.input) <- rownames(test.res)
    
    survInputs$score      <- figure.input
    survInputs$indexRow   <- 1
    survInputs$covariates <- make_covariate_df(names(figure.input), input$metagenes)
    
    att$set(85, "Generating outputs")
    
    output$update_survival_btn <- renderUI({
      actionBttn("bttn_update_surv",
                 label = tagList(icon("sync"), "Update Survival Plot"),
                 color = "warning",
                 style = "bordered",
                 size  = "sm")
    })
    
    output$Mval <- renderDT({
      datatable(test.res, options = list(pageLength = 10))
    })
    
    att$set(95, "Rendering plot")
    
    output$figure <- renderPlot({
      req(Mvals$data, Mvals$metagene)
      fi  <- setNames(Mvals$data$Risk_Value, rownames(Mvals$data))
      idx <- input$Mval_row_last_clicked
      switch(Mvals$metagene,
             "SHH"             = generate_figure_highlight_SHH(fi, idx),
             "Group3/4 (Early)"= generate_figure_highlight_G3_G4_sub(fi, idx),
             "Group3/4 (Late)" = generate_figure_highlight_G3_G4_no_sub(fi, idx))
    })
    
    att$set(100, "Complete")
    
    output$metagenechoice <- renderText(input$metagenes)
    output$sample         <- renderText(input$Mval_cell_clicked$value)
    
    showTab("tabs", "Results", select = TRUE)
    showTab("tabs", "Download")
    
    output$time <- renderText(proc.time() - ptm)
  })
  
  # ── Reset ──────────────────────────────────────────────────────────────────
  
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = TRUE)
    session$reload()
  })
  
  # ── Download ───────────────────────────────────────────────────────────────
  
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
        req(survInputs$score, survInputs$covariates)
        
        score    <- survInputs$score
        covs     <- survInputs$covariates
        indexRow <- survInputs$indexRow
        metagene <- input$metagenes
        
        # Rebuild the risk plot
        figure.input        <- score
        names(figure.input) <- names(score)
        
        risk_plot <- if (metagene == "SHH") {
          generate_figure_highlight_SHH(figure.input, indexRow)
        } else if (metagene == "Group3/4 (Early)") {
          generate_figure_highlight_G3_G4_sub(figure.input, indexRow)
        } else {
          generate_figure_highlight_G3_G4_no_sub(figure.input, indexRow)
        }
        
        # Rebuild the survival plot
        surv_plot <- render_survival_plot(
          metagene = isolate(metagene),
          score    = score,
          covs     = covs,
          indexRow = indexRow
        )
        
        # Selected sample survival value for annotation
        surv_label <- tryCatch({
          if (metagene == "SHH") {
            fit <- survfit(.shh_fit,
                           newdata = data.frame(
                             fixedMG       = as.numeric(score[indexRow]),
                             ConsensusMYCN = as.numeric(covs$mycn[indexRow])
                           ))
            val <- as.numeric(summary(fit, time = 5)$surv)
            paste0("Selected sample: ", covs$sample[indexRow],
                   " — Estimated 5-year survival: ", round(val * 100, 1), "%")
            
          } else if (metagene == "Group3/4 (Early)") {
            new_row <- data.frame(
              fixedMG      = as.numeric(score[indexRow]),
              M._versus_M. = as.numeric(covs$mets[indexRow]),
              ConsensusMYC = factor(
                paste0("ConsensusMYC=", covs$myc[indexRow]),
                levels = .g34early_fit$xlevels[["strata(ConsensusMYC)"]]
              )
            )
            fit <- survfit(.g34early_fit, newdata = new_row)
            val <- as.numeric(summary(fit, time = 5)$surv)
            paste0("Selected sample: ", covs$sample[indexRow],
                   " — Estimated 5-year survival: ", round(val * 100, 1), "%")
            
          } else {
            fit <- survfit(.g34late_fit,
                           newdata = data.frame(
                             fixedMG = as.numeric(score[indexRow])
                           ))
            val <- as.numeric(summary(fit, time = 10)$surv)
            paste0("Selected sample: ", covs$sample[indexRow],
                   " — Estimated 10-year survival: ", round(val * 100, 1), "%")
          }
        }, error = function(e) "Survival estimate unavailable")
        
        pdf(file, width = 14, height = 12)
        
        # Page 1: results table
        grid.newpage()
        gridExtra::grid.table(round(Mvals$data, 3))
        
        # Page 2: risk plot + survival plot arranged vertically
        grid.newpage()
        print(ggarrange(
          risk_plot,
          surv_plot,
          ncol  = 1,
          nrow  = 2,
          labels = c("Risk Score Distribution", "Survival Estimate")
        ))
        
        # Annotation below plots
        grid.newpage()
        grid.text(surv_label, x = 0.5, y = 0.5,
                  gp = gpar(fontsize = 14, fontface = "bold"))
        
        # Covariate table page
        if (metagene != "Group3/4 (Late)") {
          grid.newpage()
          title_grob <- grid.text("Sample Covariates", x = 0.5, y = 0.95,
                                  gp = gpar(fontsize = 14, fontface = "bold"),
                                  draw = FALSE)
          table_grob <- gridExtra::tableGrob(covs, rows = NULL)
          gridExtra::grid.arrange(title_grob, table_grob, ncol = 1,
                                  heights = c(0.1, 0.9))
        }
        
        dev.off()
      }
    }
  )
}

shinyApp(ui, server)
