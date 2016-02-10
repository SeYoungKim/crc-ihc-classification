library(shiny)
library(randomForest)
load("classifier.RData")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
   datasetInput <- reactive({
    a1=data.frame(CDX2_StainInt=input$CDX2, #levels=c("low","mod", "peak")),
               HTR2B_StainInt=input$HTR2B, #levels=c("low", "low", "mod")),
               ZEB1_StAreaFrac.norm=input$ZEB1, #factor(input$ZEB1, levels=c("low", "mod")),
               KER_StainInt=input$KERS, 
              CDX2_StAreaFrac=input$CDX2C*input$KERC/100,
              FRMD6_StainInt=input$FRMD6I, #factor(input$ZEB1, levels=c("low", "mod")),
              CDX2_StAreaFrac.norm=input$CDX2C,
              FRMD6_StAreaFrac=input$FRMD6C*input$KERC/100, 
               FRMD6_StAreaFrac.norm=input$FRMD6C, 
                             KER_StAreaFrac=input$KERC
                      )
    levels(a1$CDX2_StainInt)=list("low"="low", "mod"="mod", "peak"="high")
    levels(a1$HTR2B_StainInt)=list("low"="low", "low"="mod", "mod"="high")
   levels(a1$FRMD6_StainInt)=list("low"="low", "mod"="mod", "mod"="high")
    levels(a1$KER_StainInt)=list("low"="low", "mod"="mod", "peak"="high")
    levels(a1$ZEB1_StAreaFrac.norm)=list("low"="absent", "mod"="present")
    a1$CDX2_StAreaFrac=ifelse(a1$CDX2_StAreaFrac<5, "low", ifelse(a1$CDX2_StAreaFrac>20, "peak", "mod"))
    a1$CDX2_StAreaFrac=factor(a1$CDX2_StAreaFrac, levels=c("low", "mod", "peak"))
    a1$CDX2_StAreaFrac.norm=ifelse(a1$CDX2_StAreaFrac.norm<30, "low", ifelse(a1$CDX2_StAreaFrac.norm>50, "peak", "mod"))
    a1$CDX2_StAreaFrac.norm=factor(a1$CDX2_StAreaFrac.norm, levels=c("low", "mod", "peak"))
    
   a1$FRMD6_StAreaFrac=ifelse(a1$FRMD6_StAreaFrac<5, "low", ifelse(a1$FRMD6_StAreaFrac>40, "peak", "mod"))
    a1$FRMD6_StAreaFrac=factor(a1$FRMD6_StAreaFrac, levels=c("low", "mod", "peak"))
    a1$FRMD6_StAreaFrac.norm=ifelse(a1$FRMD6_StAreaFrac.norm<70, "low", "mod")
    a1$FRMD6_StAreaFrac.norm=factor(a1$FRMD6_StAreaFrac.norm, levels=c("low", "mod")) 
   a1$KER_StAreaFrac=ifelse(a1$KER_StAreaFrac<40, "low", "mod")
    a1$KER_StAreaFrac=factor(a1$KER_StAreaFrac, levels=c("low", "mod"))
    a1
    })
  
 #
  
  output$summary <- renderPrint({
    dataset <- datasetInput()
    dataset
  })
  
  
  output$Classification <-renderUI({
    ClassRslt=predict(rf.test, datasetInput(), type="prob")
    str1 <- paste("\n Probability of epithelial-like subtype: ", ClassRslt[1])
    str2 <- paste("\n Probability of mesenchymal-like subtype:",  ClassRslt[2])
    HTML(paste(str1, str2, sep = '<br/>'))
  })
  
 
  output$view <- renderTable({
    head(datasetInput())
  })
  
  output$Nplot <-renderPlot({
    ClassRslt=predict(rf.test, datasetInput(), type="prob")
    barplot(t(ClassRslt), horiz=T, beside=F, col=c("lightblue", "darkblue"), names.arg=" ",
            legend.text=c("epithelial-like", "mesenchymal-like"), xlab="prediction probability")
    }
    )
})
