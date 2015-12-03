function [BMap, BSum] = BrownMap(imHED, X)
% BROWNMAP: Selects brown regions using Otsu thresholding
%
%   [BMap, BSum] = BrownMap(imHED) takes a Colour Deconvoluted TMA image (imHED) 
%   processed by ColourDeconvolve to generate a Otsu thresholded Map of DAB
%   regions in the image (Map) and a numerical Output (AvInt) of the 
%   Average grayscale Brown Intensity.
%
%   [BMap, BSum] = BrownMap(imHED, X) performs this operation on the Xth
%   dimension of this image.
%
%   Version 2. Anne Trinh
    tic
    % check inputs
    N=nargin;
    error(nargchk(1,2,N));
    if ~exist('X', 'var'), X=3; end
    assert(ndims(imHED)>=X, 'Number of channels of input image does not match the specified image dimension');
    
    % process the DAB channel and select the pixels with intensity<255 
    % (ie. pixels which constitute the TMA)
    imDAB=imHED(:,:,X);
%     DABvec=imDAB(:);
%     DABvec=im2double(DABvec)*255;
%     x1=find(DABvec>0);
%     DABvec=DABvec(x1);
    BMap=zeros(size(imDAB));

    % Find an optimal threshold value using the otsu method on selected
    % pixels
    t1=graythresh(imDAB);
    if t1<0.35, t1=0.35; end
    
    BMap=im2bw(imDAB, t1);
    %BMap=~BMap;
    x2=find(BMap==1);
  %  PropS=length(x2)/length(x1);
    MInt=mean(imDAB(x2));
    BSum=struct('StainArea', length(x2),'MeanIntensity', MInt, 'Threshold', t1); % , 'QuantAllredScore', 5*PropS+3*(255-MInt)/255);
%   BSum=[length(x2), MInt];
   fprintf(sprintf('\n'))
toc
  
