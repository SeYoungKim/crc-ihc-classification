##Collate the information from Image Analysis into a single R object
# Here we use the AMC set as an example:

source("classifier_functions.R")

#########################
## 1. Prepare the data
########################
# Firstly run a script to ensure the data and the maps are aligned with the same patient:
# Requires a folder containing a .csv document containing the stain files of interest, as well 
# as a folder containing the mapping information. The "mapping" folder must contain:
# 1. a file containing "_map.csv" extension. This is the TMA map. May contain several files if several slides are present
# 2. "_slide_no.csv" file containing unique identifiers to determine which 
ReadFromCSV("AMC_data/", "AMC_mapping/")

# The output of this is an .RData file. Let's load this:
AMCdata=LoadFiles("AMC_data")

# Next, we need to load information about the patient, including MSI status and class if this is the predictor class.
AMCclin=read.csv("AMCclinical.csv", header=T)

#Next, we want to remove patients with missing cores and MSI positive (or unknown) patients to produce a two-class classification system (between CCS1 and CCS3).
#Inputs: *the data, *the clinical information, *Train or Test data *TMA area estimate: slides containing values beyond this are discarded. It is expected
# that a full core would occupy 30-60% of an image
AMCdata.train=PreProcess(AMCdata, AMCclin[ , c("ID", "MSI2", "Class")],
                           'Train', c(0.15, 0.6)) 

# This produces an object with a $data list (containing image information) and a $PatID list (containing Patient ID, MSI, CCS type, Image No etc)

#########################
## 2. Create Classifier
########################
## Here we train the AMC classifier
set.seed(1119)
## Define the true labels based on the predicted class
TrueLab=factor(AMCdata.train$PatID$predClass)
## Prepare the training data by scaling, and removing information about the image ("ImageSize" is used only for QC)
AMC_raw=MedScale(AMCdata.train$data[, -grep("Im", colnames(AMCdata.train$data))])
tdat <- data.frame(AMC_raw, CCS=factor(TrueLab))
## Train the classifier!
rf.TMA <- randomForest(CCS ~ ., data=tdat, importance=TRUE, ntree=1000, proximity=T)
rf.TMA

#########################
## 3. Predict Output
########################
# use the classifier to predict a patient cohort:
output = predict(rf.TMA, AMC_raw, type="prob")
# Collapse to a single patient score:
PatClasses = Compress.Scores(output, AMCdata.train$PatID$Pat, 2, 0.6)

#########################
## 4. Example with CairoSlides
########################

## Run the above loading steps to classify this subset of CAIRO patients
ReadFromCSV("CAIRO_data", "CAIRO_mapping")
Cai1data=LoadFiles("CAIRO_data/")
## Input clinical data
Cai1Clin=read.csv("Cairo1_clinical.csv")
## preprocess data
Cai1data.edit=PreProcess(Cai1data, Cai1Clin[ , c("ID", "MSI")], 'Test', c(0.1, 0.8))
Cai1_raw=Cai1data.edit$data[, -grep("Im", colnames(Cai1data.edit$data))]
                          
## Prepare and predict!
Cai1_range=MedScale(Cai1_raw)
# Make the prediction
C1_R<- predict(rf.TMA, Cai1_range, type="prob")
# Compress the Cores
C1_RPc=Compress.Scores(C1_R, Cai1data.edit$PatID$Pat, 2, 0.6)
