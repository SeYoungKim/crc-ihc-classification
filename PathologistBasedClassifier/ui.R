library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Colorectal Cancer Subtypes by IHC"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("CDX2", "CDX2 Stain Intensity", 
                  choices = c("low","mod", "high")),
      selectInput("HTR2B", "HTR2B Stain Intensity", 
                  choices = c("low", "mod", "high")),
    selectInput("FRMD6I", "FRMD6 Stain Intensity", 
                choices = c("low", "mod", "high")),
      selectInput("ZEB1", "ZEB1 Epithelial content", 
                  choices = c("absent", "present")),
      selectInput("KERS", "KER Stain Intensity", 
                  choices = c("low","mod", "high")),
    sliderInput("CDX2C", "CDX2 content", 0, 100,
                0, 5),
    sliderInput("FRMD6C", "FRMD6 content", 0, 100,
                50, 5),
    sliderInput("KERC", "Keratin content", 0, 100,
                50, 5),
      submitButton()
    ),
  
    
      
   mainPanel(
   tabsetPanel(
     tabPanel("Compute Classification", h3("Classification result:"),
              h4(htmlOutput("Classification")),
              plotOutput("Nplot")),
     tabPanel("CDX2", h4("Examples of CDX2 scoring:"),
              img(src="cdx2_content.jpg"),
              h5("Examples of epithelial specific stain content: 1. 25% 2. 55% 3. 95%"),
              img(src="cdx2_intensity.jpg"),
              h5("Examples of 1. Low (Left), 2. Medium (Middle), 3. High (Right)  intensity staining ")
               ), 
     tabPanel("FRMD6", h4("Examples of FRMD6 scoring:"),
              img(src="frmd6_content.jpg"),
              h5("Examples of epithelial specific stain content: 1. 10% 2. 60% 3. 90%"),
              img(src="frmd6_intensity.jpg"),
              h5("Examples of 1. Low (Left), 2. Medium (Middle), 3. High (Right)  intensity staining ")),
     tabPanel("HTR2B", h4("Examples of HTR2B scoring:"),
              img(src="htr2b_intensity.jpg"),
              h5("Examples of 1. Low (Left), 2. Medium (Middle), 3. High (Right)  intensity staining ")),
     tabPanel("ZEB1", h4("Examples of ZEB1 scoring:"),
              img(src="zeb1_intensity.jpg"),
              h5("Examples of 1. Absent (Left) and 2. Epithelial specific scoring (Right")),
     tabPanel("KER", h4("Examples of KER scoring:"),
              img(src="ker_content.jpg"),
              h5("Examples of keratin content within a core: 1. 15% 2. 50% 3. 85%"),
              img(src="ker_intensity.jpg"),
              h5("Examples of 1. Low (Left), 2. Medium (Middle), 3. High (Right)  intensity staining ")
             )
   ))
  )
))

