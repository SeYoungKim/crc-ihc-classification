%% Run image analysis Script
% 	This script contains the methods run on page 3-4 of the supplementary information

%Set the working directory and name the file of interest
fprintf('Load the image of interest')
filename ='R12C3' ;

% Read in the image
ImA = imread(sprintf('%s.tif', filename));
%Show the image size 
size(ImA)

% Step1 . Find the main TMA Area:
[imA2, TMAStat, TMAOut]=TMAOutline(ImA); 

% Information on the core is save in 'TMAStat'
TMAStat

% Step 2. Colour Deconvolve the image 
imHD=ColourDeconvolve(imA2);

% Step 3. Detect the Brown Region
[imBrown, BStat]=BrownMap(imHD);
BStat

% Step 4. Save the output
% write the adjusted image
imwrite(imA2, sprintf('~/%sTMA.jpg', filename), 'jpeg');
% Highlight the main TMA area
imwrite(SegArea(imA2, TMAOut), sprintf('~/%sOutline.jpg', filename), 'jpeg');
% Colour Deconvolution
imwrite(imHD, sprintf('~/%sCD.jpg', filename), 'jpeg');
% Highlight the "brown" region
imwrite(SegArea(imA2, imBrown), sprintf('~/%sBrown.jpg', filename), 'jpeg');


% write the following for classification in R:
StatsOut={TMAStat.TotalArea, TMAStat.PercentImage, BStat.StainArea, BStat.MeanIntensity, BStat.StainArea/TMAStat.TotalArea};
StatsOut