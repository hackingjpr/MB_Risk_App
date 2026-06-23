source("./MB_app_source_functionsJun26.R")


ui <- shiny::fluidPage(
  
  #Loading bar
  useAttendant(),
  #Dashboard looks etc
  dashboardPage(
    skin = "blue",
    dashboardHeader(title= "", dropdownMenuOutput("messageMenu")), 
                dashboardSidebar(
                  #logo
                  # (img(src='Free_Sample_By_Wix (5).jpg', align = "center")),
                  #(img(src='Free_Sample_By_Wix%20(5).jpg', align = "center")),
                  (img(src='MB_risk.jpg', align = "center")),
                  #uploading idat
                  h2("Step 1"),
                  fileInput(
                    "idatFile",
                    "Upload idat file:",
                    multiple = TRUE,
                    accept = ".idat"),
                  #Step 2
                  h2("Step 2"),
                  selectInput("metagenes", "Select Metagene Set", c("SHH", "G3_G4_no_sub", "G3_G4_sub")),
                  #set it away
                  h2("Step 3"),

                  
                  actionBttn(
                    inputId = "bttn1",
                    label = "Generate Risk Values",
                    color = "primary",
                    style = "jelly"
                  ),
                  
                  
                  #Loading Bar
                  attendantBar("progress-bar"),
                  #Reset button
                  column(width = 12,
                         "Want to upload a different data set?"),
                  column(
                    width = 12,
                    actionBttn(
                      inputId = "bttn2",
                      label = "Reset",
                      color = "danger",
                      style = "fill",
                      size = "sm"
                    ),
                    align = "center"
                  ),
                  
                  hr(),
                  hr(),
                  hr(),
                  hr(),

                  tags$footer(
                    hr(),
                    HTML(
                    paste(
                      "Author: James Hacking,", "<br/>",  
                      "Date Created: 12-06-2026,", "<br/>", 
                      "Copyright (c) James Hacking, 2026,", "<br/>",
                      h5(
                        a("Email: james.hacking@ncl.ac.uk", href = "mailto:james.hacking@ncl.ac.uk")
                      ),
                      h4(
                        "Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic or clinical use."
                      )
                      
                    )))
                ), dashboardBody((
                  tabsetPanel( id = "tabs",
                               #Information Tab
                               tabPanel(
                                 "Info",
                                 (fluidRow(
                                   includeMarkdown("./introduction.md")
                                  )
                                 )
                               ),
                               #Tutorial Tab
                               tabPanel("Tutorial",
                                        includeMarkdown("./Tutorial/tutorial.md")
                                          ),
                               #Results Tab
                               tabPanel("Results Table",
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Risk Values",
                                            status = "info",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            DTOutput('Mval')
                                          ))),
                                        (fluidRow(
                                          box(
                                            title = "Risk plot", 
                                            status = "warning", 
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("figure")),
                                          box(
                                            title = "Selections",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            h3("Metagene Set:"),
                                            textOutput("metagenechoice"),
                                            
                                            h3("Sample Selected:"),
                                            textOutput("sample"),
                                            
                                            h3("Patient's Risk Percentile:"),
                                            textOutput("percentages"),
                                            
                                            h3("Patient's Survival Percentile:"),
                                            textOutput("survivalPercentage"),
                                            
                                            h4("Disclaimer : This app is designed exclusively for research purposes and is strictly not for diagnostic use.")
                                          ))
                                        ),
                                        (fluidRow(
                                          box(
                                            #width = 12,
                                            title = "Survival Curve With No Variables Considered",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("figureSurvival")
                                          ),
                                          box(
                                            #width = 12,
                                            title = "Survival Curve With Age Considered",
                                            status = "warning",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            plotOutput("figureSurvivalAge")
                                          )
                                        )),
                                        (fluidRow(
                                          box(
                                            width = 12,
                                            title = "Summary",
                                            status = "success",
                                            solidHeader = TRUE,
                                            collapsible = TRUE,
                                            #includeMarkdown("./AppExtraFiles/ResultsSummary.md")
                                          )
                                        ))
                                        ),
                               
                               #Download Tab
                               tabPanel("Download",
                                        textInput("filename", "Please insert desired filename", "M-values"),
                                        radioButtons(inputId = "download", label = "Select file type", choices = c("csv", "pdf")),
                                        downloadButton("down", "Download the results"),
                               )
                  )
                )),
                #Font Selection            
                tags$head(tags$style(HTML('* {font-family: "Courier New"};')))),
  
)

server <- function(session, input, output) {
  # Allow larger files to be uploaded
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  #Hide tabs until graph made
  hideTab(inputId = "tabs", target = "Results Table")
  hideTab(inputId = "tabs", target = "Download")
  
  #Loading bar
  att <- Attendant$new("progress-bar")
  unlink("./temp/", recursive = T)
  if(!dir.exists("./temp")){
    dir.create("./temp")
  }
  #fix the overwriting issue
  fs::path("./temp/",createRandString()) -> tempDIR
  dir.create(tempDIR)
  
  Mvals <- reactiveValues(data = NULL)
  observeEvent(input$bttn1, {
    
    req(input$idatFile)
    
    att$set(0, text = "Starting...")
    
    # STEP 1 — Copy files
    input$idatFile$datapath -> in.files
    paste0(tempDIR, "/", input$idatFile$name) -> out.files
    
    file.copy(in.files, out.files)
    att$set(15, text = "Files uploaded")
    
    # STEP 2 — Get basenames
    temp.base <- get_basenames(tempDIR)
    att$set(25, text = "Preparing data")
    
    # STEP 3 — Process IDAT files (HEAVY STEP)
    ptm <- proc.time()
    temp.processed <- process_idats(temp.base)
    att$set(60, text = "Processing methylation data")
    
    # STEP 4 — Select metagene
    meta <- if (input$metagenes == "SHH") {
      SHH
    } else if (input$metagenes == "G3_G4_sub") {
      G3_G4_sub
    } else {
      G3_G4_no_sub
    }
    
    att$set(70, text = "Applying metagene model")
    
    # STEP 5 — Calculate risk
    test.res <- extract.metagene(
      as.character(meta[[1]]$genes),
      as.numeric(meta[[1]]$weights),
      beta2m(temp.processed$betas),
      as.numeric(meta[[2]])
    )
    
    test.res <- round(test.res, 3)
    att$set(85, text = "Generating outputs")
    
    # STEP 6 — Render outputs
    output$Mval <- renderDT({
      test.res
    })
    
    figure.input <- test.res$Risk_Value
    names(figure.input) <- rownames(test.res)
    
    att$set(95, text = "Rendering plots")
    
    output$figure <- renderPlot({
      if (input$metagenes == "SHH") {
        generate_figure_highlight_SHH(figure.input, input$Mval_row_last_clicked)
      } else if (input$metagenes == "G3_G4_sub") {
        generate_figure_highlight_G3_G4_sub(figure.input, input$Mval_row_last_clicked)
      } else {
        generate_figure_highlight_G3_G4_no_sub(figure.input, input$Mval_row_last_clicked)
      }
    })
    
    # FINAL STEP
    att$set(100, text = "Complete")
    
    showTab(inputId = "tabs", target = "Results Table", select = TRUE)
    showTab(inputId = "tabs", target = "Download")
    
    output$time <- renderText({ proc.time() - ptm })
    
  })
  
  #Reset session and delete tempDIR
  observeEvent(input$bttn2, {
    unlink(tempDIR, recursive = T)
    session$reload()
    return()
    print("session reload not working")
  })
  renderText(output$metagenes <- input$metagenes)
  # output$down <- downloadHandler(
  #   filename = function() {
  #     paste(input$filename,Sys.time(), ".csv", sep="_")
  #   },
  #   content ={ 
  #     function(file) {
  #       write.csv(Mvals(), file, row.names = TRUE)
  #     }
  #     # if (input$download == "csv")
  #     #   write.csv(output$Mval)
  #     # else
  #     #   pdf(output$Mval,
  #     #       width = 14
  #     #   )
  #     # dev.off()
  #   })
  output$down <- downloadHandler(
    filename = function() {
      paste(input$filename, Sys.time(), ifelse(input$download == "csv", ".csv", ".pdf"), sep = "_")
    },
    content = function(file) {
      if (input$download == "csv") {
        write.csv(Mvals$data, file, row.names = TRUE)
      } else if (input$download == "pdf") {
        pdf(file, width = 14, height = 8)
        
        Mvals$data <- round(Mvals$data, digits =3)
        grid.table(Mvals$data)
        # 
        # message("input")
        # input$idatFile$datapath -> in.files
        # message(in.files)
        # paste0(tempDIR,"/", input$idatFile$name) -> out.files
        # message(out.files)
        # 
        # copied <- file.copy(in.files, out.files)
        # message("Copied")
        # 
        # temp.base <- get_basenames(tempDIR)
        # message("get_basenames done")
        # cat("Timing start\n")
        # ptm <- proc.time()
        # 
        # temp.processed <- process_idats(temp.base)
        # on.exit({
        #   att$done()
        # })
        # message("process_idats done")
        # # test.res <-(  
        # #   
        # # if (input$metagenes == "MRT (ATRT & ECRT)") {
        # #     ALL -> meta
        # #   } %>%
        # # if (input$metagenes == "ATRT") {
        # #   ATRT -> meta
        # # } %>%
        # # if (input$metagenes == "ECRT") {
        # #   ECRT -> meta
        # # } %>%
        # # #this needs to be reactive
        # #  #output$test.res <- reactive({extract.metagene(
        # #  extract.metagene(
        # #   as.character(meta[[1]]$genes),
        # #   as.numeric(meta[[1]]$weights),
        # #   beta2m(temp.processed$betas),
        # #   as.numeric(meta[[2]])
        # # ))
        # #)}))
        # meta <- reactive({
        #   if (input$metagenes == "MRT (ATRT & ECRT)") {
        #     return(ALL)
        #   } else if (input$metagenes == "ATRT") {
        #     return(ATRT)
        #   } else if (input$metagenes == "ECRT") {
        #     return(ECRT)
        #   } else {
        #     # Handle other cases if needed
        #     return(NULL)
        #   }
        # })
        # 
        # # test.res <- reactive({
        # test.res <- if (!is.null(meta())) {
        #   extract.metagene(
        #     as.character(meta()[[1]]$genes),
        #     as.numeric(meta()[[1]]$weights),
        #     beta2m(temp.processed$betas),
        #     as.numeric(meta()[[2]])
        #   )
        # }
        # # })
        # message("test.res done")
        # message(test.res)
        # # as.numeric(test.res)
        # round(test.res, digits = 3) -> test.res
        # output$resultText <- renderText({
        #   # Display the result in the UI
        #   req(test.res())  # Ensure that test.res() is not NULL
        #   result <- round(test.res(), digits = 3)
        #   paste("Result: ", result)
        # })
        # 
        # message("round test.res done")
        # output$Mval <- renderDT (({test.res
        # }),
        # options = list(
        #   pageLength = 10, 
        #   processing=FALSE),
        # selection = "single"
        # )
        # message("render dt done")
        # print(test.res)
        # figure.input <- test.res$Risk_Value
        # names(figure.input) <- rownames(test.res)
        # print(figure.input)
        # 
        # att$done(text = "Complete")
        # 
        # 
        # output$time <- renderText({proc.time() - ptm})
        # rowSelect <- reactive({input$Mval_rows_selected})
        # 
        #  x <-
        #   grid.draw(
        #     
        #     # generate_figure_highlight_mrt(
        #     #   figure.input
        #     #   ,1)
        #     
        #     generate_figure_highlight_mrt(
        #       figure.input
        #       ,input$Mval_row_last_clicked)
        #     # 
        #     # if (input$metagenes == "MRT (ATRT & ECRT)") {
        #     #   generate_figure_highlight_mrt(
        #     #     figure.input
        #     #     ,input$Mval_row_last_clicked)
        #     # },
        #     # elseif (input$metagenes == "ATRT") {
        #     #   generate_figure_highlight_atrt(
        #     #     figure.input
        #     #     ,input$Mval_row_last_clicked)
        #     # },
        #     # elseif (input$metagenes == "ECRT") {
        #     #   generate_figure_highlight_ecrt(
        #     #     figure.input
        #     #     ,input$Mval_row_last_clicked)
        #     # }
        #     
        #     # https://rstudio.github.io/DT/shiny.html
        #     
        #     
        #   )
        #  x
        # 
        # # Print DataTable
        # # print(DT::datatable(Mvals$data))
        # # 
        # # # Print ggplot figure
        # # print(generate_figure_highlight_mrt(Mvals$data, 1))
        # 
        
        dev.off()
      }
    }
  )
  
  
  
  
  
}

shinyApp(server = server, ui = ui)
