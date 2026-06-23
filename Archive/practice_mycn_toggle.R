ui <- fluidPage(
  checkboxInput("mycn_toggle", "Set MYCN = 1", value = FALSE),
  actionButton("update_plot", "Redraw Plot"),
  plotOutput("myPlot"),
  selectInput("sample_select", "Choose sample", choices = paste0("S", 1:10))
)

server <- function(input, output, session) {
  
  # Original dataframe (never changes)
  base_df <- reactive({
    data.frame(
      sample = paste0("S", 1:10),
      MYCN = 0
    )
  })
  
  # Reactive dataframe that updates based on checkbox
  modified_df <- reactive({
    df <- base_df()
    
    if (input$mycn_toggle) {
      df$MYCN[df$sample == input$sample_select] <- 1
    }
    
    return(df)
  })
  
  # Only redraw plot when button is clicked
  plot_data <- eventReactive(input$update_plot, {
    modified_df()
  })
  
  output$myPlot <- renderPlot({
    df <- plot_data()
    
    barplot(df$MYCN, names.arg = df$sample)
  })
}


shinyApp(ui, server)
