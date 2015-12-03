library(shiny)
library(randomForest)
load("classifier.RData")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  
  ## Note: should probably group low and moderate staining samples together?
  datasetInput <- reactive({
    a1=data.frame(CDX2_StainInt=input$CDX2, #levels=c("low","mod", "peak")),
               HTR2B_StainInt=input$HTR2B, #levels=c("low", "low", "mod")),
               ZEB1_StAreaFrac.norm=input$ZEB1, #factor(input$ZEB1, levels=c("low", "mod")),
               KER_StainInt=input$KERS, #factor(input$KERS, levels=c("low","mod", "peak")),
              CDX2_StAreaFrac=input$CDX2C*input$KERC/100,
              #factor(cut(input$CDX2C*input$KERC/100, c(-1, 10, 30, 105) , levels=c("low", "mod", "peak"))),
              CDX2_StAreaFrac.norm=input$CDX2C,
              #factor(cut(input$CDX2C, c(-1, 50, 70, 101), levels=c("low", "mod", "peak"))),
              FRMD6_StAreaFrac=input$FRMD6C*input$KERC/100, #factor(cut(input$FRMD6C*input$KERC/100, c(-1, 45, 101), levels=c("low", "mod"))),
              FRMD6_StAreaFrac.norm=input$FRMD6C, #factor(cut(input$FRMD6C, c(-1, 70, 101), levels=c("low", "mod")))
              KER_StAreaFrac=input$KERC #, levels=c("low", "mod"))
                      )
   # browser()
    levels(a1$CDX2_StainInt)=list("low"="low", "mod"="mod", "peak"="high")
   # browser()
    levels(a1$HTR2B_StainInt)=list("low"="low", "low"="mod", "mod"="high")
    levels(a1$KER_StainInt)=list("low"="low", "mod"="mod", "peak"="high")
    levels(a1$ZEB1_StAreaFrac.norm)=list("low"="absent", "mod"="present")
    a1$CDX2_StAreaFrac=ifelse(a1$CDX2_StAreaFrac<5, "low", ifelse(a1$CDX2_StAreaFrac>20, "peak", "mod"))
    a1$CDX2_StAreaFrac=factor(a1$CDX2_StAreaFrac, levels=c("low", "mod", "peak"))
    a1$CDX2_StAreaFrac.norm=ifelse(a1$CDX2_StAreaFrac.norm<30, "low", ifelse(a1$CDX2_StAreaFrac.norm>50, "peak", "mod"))
    a1$CDX2_StAreaFrac.norm=factor(a1$CDX2_StAreaFrac.norm, levels=c("low", "mod", "peak"))
    
    a1$FRMD6_StAreaFrac=ifelse(a1$FRMD6_StAreaFrac<45, "low", "mod")
    a1$FRMD6_StAreaFrac=factor(a1$FRMD6_StAreaFrac, levels=c("low", "mod"))
    a1$FRMD6_StAreaFrac.norm=ifelse(a1$FRMD6_StAreaFrac.norm<70, "low", "mod")
    a1$FRMD6_StAreaFrac.norm=factor(a1$FRMD6_StAreaFrac.norm, levels=c("low", "mod"))
    a1$KER_StAreaFrac=ifelse(a1$KER_StAreaFrac<40, "low", "mod")
    #browser()
    a1$KER_StAreaFrac=factor(a1$KER_StAreaFrac, levels=c("low", "mod"))
    a1
    })
  
 #
  
  output$summary <- renderPrint({
    dataset <- datasetInput()
    dataset
  })
  
  
  output$Classification <-renderUI({
  #  browser()
    ClassRslt=predict(rf.test, datasetInput(), type="prob")
    #output$CCS2=ClassRslt[3]
    #ClassRslt
    #sprintf(paste("\n Probability of epithelial-like subtype: ", ClassRslt[1],
  #          "\n Probability of mesenchymal-like subtype:",  ClassRslt[2]))
    str1 <- paste("\n Probability of epithelial-like subtype: ", ClassRslt[1])
    str2 <- paste("\n Probability of mesenchymal-like subtype:",  ClassRslt[2])
    HTML(paste(str1, str2, sep = '<br/>'))
  })
  
 
  output$view <- renderTable({
    head(datasetInput())
  })
  
})
