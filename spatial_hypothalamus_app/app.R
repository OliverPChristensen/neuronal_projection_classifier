# Load the shiny package
library(shiny)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Simple Shiny App"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      # Input: Slider for the number of observations
      sliderInput("obs", 
                  "Number of observations:", 
                  min = 1, 
                  max = 1000, 
                  value = 500)
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      # Output: Plot of randomly generated data
      plotOutput("plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Generate a plot of random data
  output$plot <- renderPlot({
    data <- rnorm(input$obs)
    hist(data, main = "Histogram of Random Data")
  })
}

# Run the application
shinyApp(ui = ui, server = server)