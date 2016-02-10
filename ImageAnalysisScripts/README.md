# crc-ihc-classification

This folder contains the scripts and instructions to run the automated classifier.

Prior to this, the following are required:
* TMA cores segmented into a series of image cores
* A TMA map containing patient information
* Installation of MATLAB (v7.14 or later) and R

IMAGE ANALYSIS:
An example of the image segmentation pipeline is scripted in TMAPipeline/scriptRun.m for one given core.
This can be adapted to loop through multiple images.
A csv file containing a matrix of 5 features per image is required for classification.

CLASSIFICATION:
An example of the classification approach using the AMC training set is in classifier/RunClassifier.R
This classifier is to be applied to the CAIRO dataset.
For each cohort, three items are required: 
(1) 5 x csv files with results from MATLAB analysis (see AMC_data for examples)
(2) maps. At least two files required: 1. TMA maps 2. a Table with identifiers: which slide corresponds to which stain? (see AMC_mapping/C1_slide_no.csv)
(3) Clinical Information. Patient IDs with MSI status to determine CMS1/MSI+ patients.

