library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Colorectal Cancer Subtypes by IHC"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    mainPanel(
      selectInput("CDX2", "CDX2 Stain Intensity", 
                  choices = c("low","mod", "high")),
     # selectInput("CDX2A", "CDX2 Epithelial content", 
    #              choices = c("low", "mod", "peak")),
     # selectInput("FRMD6", "FRMD6 Epithelial content", 
    #              choices = c("low", "mod", "high")),
      selectInput("HTR2B", "HTR2B Stain Intensity", 
                  choices = c("low", "mod", "high")),
      selectInput("ZEB1", "ZEB1 Epithelial content", 
                  choices = c("absent", "present")),
      selectInput("KERS", "KER Stain Intensity", 
                  choices = c("low","mod", "peak")),
    sliderInput("CDX2C", "CDX2 content", 0, 100,
                50, 5),
    sliderInput("FRMD6C", "FRMD6 content", 0, 100,
                50, 5),
    sliderInput("KERC", "Keratin content", 0, 100,
                50, 5),
    #  selectInput("KERA", "KER Stain Content", 
    #              choices = c("low", "mod")),
    #  helpText("Note: while the data view will show only the specified",
    #           "number of observations, the summary will still be based",
    #           "on the full dataset."),
      submitButton()
    ),
  
    
      
    # Show a plot of the generated distribution
    sidebarPanel(
    # h3(textOutput("caption", container = span)),
   #  print("Input data:"),
    # verbatimTextOutput("summary"), 
   #  tableOutput("view"),
    h3("Classification result:"),
    h4(htmlOutput("Classification"))
    #h4(sprintf("Probability of mesenchymal-like phenotype: %s")),
   # tableOutput("Classification")
    )
  )
))

