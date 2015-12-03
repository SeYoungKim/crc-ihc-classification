# crc-ihc-classification

This folder contains the pathologist friendly classifier developed. 

To run in RStudio:
>> install.packages("shiny")
>> install.packages("randomForest")
>> library(shiny)
>> runApp("PathologistBasedClassifier")

Required inputs:
* CDX2 stain intensity (absent to low / moderate / high)
* CDX2 epithelial specific content as a rough percentage
* FRMD6 epithelial specific content as a rough percentage
* ZEB1 epithelial content as a binary input
* HTR2B stain intensity (absent to low / moderate / high)
* KERATIN stain intensity (absent to low / moderate / high)
* KERATIN content as a rough percentage

The output is a probability from random forest in classifying your sample into one of two classes

