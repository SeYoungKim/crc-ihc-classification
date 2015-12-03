function[imOut Stat SSeg] = TMAOutline(Samp, varargin)
% TMAOutline: Outline the TMA and report segment statistics
%
%     [imOut, Stats, SSeg] = TMAOutline(Samp) determines the main TMA region in the 
%     image Samp, which is returned in SSeg with background intensities set to 255. 
%     Statistics are returned in Stat, including 
%     the following properties:
%           Area, Percentage Area of the Image, Threshold values used,
%           Major and Minor Axis lengths, Mean Intensity, Solidity and the
%           Centroid. 
%     Note that the image may not contain a single continuous TMA section.
%     In this case, the area of each segment is given.  
%
%     The program detects the main image area by using an entropy filter, 
%     with thresholding set at double the modal entropy value for each image.
%
%     [imOut, Stats] = TMAOutline(Samp, 'threshVal', X) allows the user to designate 
%     a thresholding value, X, scaled between [0,1]. A higher value of 
%     threshVal results in a finer segmentation. The automatically detected 
%     threshold is usually in the range of 0.3-0.4 to allow detection of
%     the main TMA regardless of the staining intensity of the DAB.
%
%     In the case where there are multiple disconnected segments which
%     constitute the main TMA area, these are selected using the following 
%     criteria:
%           - solidity > 0.5 
%           - area greater than a size threshold, sizeVal, a fraction of
%           the total estimated TMA area.
%
%     [imOut, Stats] = TMAOutline(Samp, 'sizeVal', Y) returns an image of the 
%     selected area where segments are greater than Y*estimated TMA area.
%     sizeVal must be a value in (0, 1)
%
%     Example:
%     --------
%       imageA = imread('TMA_Sample_2.tif');
%       [imOut, Stats, SSeg] = TMAOutline(imageA);
%       imview(SegArea(imageA, SSeg));
%
%     Note: dependent on MATLAB Image Processing Toolbox and SelectedArea function
%     Version 2. Dated 15-6-2013 Anne Trinh

% check inputs
    tic
    Nopt=nargin-1;
    error(nargchk(1,7,Nopt+1), 'Incorrect number of input variables');
    assert(ndims(Samp)>1, 'Insert an image.');
    assert(round(Nopt/2)*2==Nopt, 'Error. A variable input argument is undefined');
    
% assign values to threshVal and sizeVal    
% default values
    threshVal=0.4;
    sizeVal=0.1;
    userVal=0;

if Nopt>=2
    for n=1:Nopt/2 
        switch varargin{2*n-1}
            case 'entVal'
                threshVal=varargin{2*n};
                assert(threshVal>=0 & threshVal<=1, '''entVal'' must be a value in the range [0,1].');
                userVal=1;
            case 'sizeVal'
                sizeVal=varargin{2*n};
                assert(sizeVal>=0 & sizeVal<=1, '''sizeVal'' must be a value in the range [0,1].'); 
            case 'intRange'
                mrange=varargin{2*n};
                assert(length(mrange)==2, '''intRange'' must have a min and max intensity input'); 
                assert(mrange(1)>=0 & mrange(1)<1 & mrange(2)>0 & mrange(2)<=1, 'Intensities must be in [0,1]'); 
        end
    end
end
    
% convert image to gray and scale between [0, 255]
    if ndims(Samp)==3
        OUT=rgb2gray(Samp);
    else 
        OUT=im2uint8(Samp);
    end
    
    f2=imcomplement(im2double(OUT));
    
    if ~exist('mrange', 'var')
        mrange=im2double([min(OUT(:)), max(OUT(:))]);
    end
    OUT=imadjust(OUT, mrange, []);

% find edges using entropy filter  
    seF = strel('disk', 2);
    F3=entropyfilt(OUT);
    F3=F3/max(F3(:))/2+f2/2;
    
% threshold depending on whether there is a user input or not
   
    switch userVal
        case 0
%             [X, C]=imhist(F3/max(max(F3)));
%             nX=length(X);
%             Y=find(X==max(X(2:nX-1))); 
%             Y2=max(max(F3))*2*C(Y+1);
%             if C(Y)>0.5
%                 Y2=max(max(F3))*threshVal;
%             end
            Y2=graythresh(F3);
        case 1
            Y2=threshVal;
    end
        
    BW=zeros(size(F3));
    a1=find(F3>Y2); %|OUT<max(OUT(:)).*0.9);
    BW(a1)=1;
    T1=Y2;

% dilate the edges to close gaps and remove small particles (<1000 x TMA
% size)
    est=sum(BW(:)); 
    B2=bwareaopen(BW, round(est/1000));
    BW=imdilate(B2, seF);                      
   % Imdiff=~BW;    
     
    
    %B2=bwlabel(B2);
% remove solid large 'holes'
%     A2=regionprops(B2, 'Area', 'Solidity');
%     ID=find([A2.Solidity]>0.8);    
%     [A3]=SelectedArea(ID, B2);
%     BigHoles=logical(B2)-logical(A3);
%     BigHoles=bwareaopen(BigHoles, round(est/15000));
%     SegIm=BW2-logical(B2)+BigHoles;
    BW2=bwlabel(BW);%(SegIm);
% Determine main TMA region. Must have an area >10% of the estimated TMA size.
    A2=regionprops(BW2, OUT, 'Area', 'Solidity', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'MeanIntensity');
    idx=find([A2.Area]>round(est*sizeVal));
    SSeg=SelectedArea(idx, BW2, 1, Samp);
% Non-TMA regions are erased (background white) & stretch the range
    
if ndims(Samp)==3

    imOutR=Samp(:,:,1);
    imOutG=Samp(:,:,2);
    imOutB=Samp(:,:,3);

   % if max(mrange)~=1
        imOutR=imadjust(imOutR, mrange,[0.02, 0.98]);
        imOutB=imadjust(imOutB, mrange,[0.02, 0.98]);
        imOutG=imadjust(imOutG, mrange,[0.02, 0.98]);
   % end
    imOutR(~SSeg)=255;
    imOutG(~SSeg)=255;
    imOutB(~SSeg)=255;
    imOut=cat(3, imOutR, imOutG, imOutB);
    
else
    imOut=imadjust(Samp);
    imOut(~SSeg)=255;
end
    
 
% Prepare statistics
        Stats=A2(idx);
        impix=size(BW,1)*size(BW,2);
        Stat.TotalArea=sum([Stats.Area]);
        Stat.PercentImage=round(sum([Stats.Area])*100/impix);
        Stat.AbsoluteEntropyThreshold=T1;
        Stat.RelativeEntropyThreshold=T1/max(max(F3));
        Stat.Area=[Stats.Area];
        Stat.Solidity=[Stats.Solidity];
        Stat.MajorAxisLength=[Stats.MajorAxisLength];
        Stat.MinorAxisLength=[Stats.MinorAxisLength];
        Stat.MeanIntensity=[Stats.MeanIntensity];
        Stat.Centroid=[Stats.Centroid];
        fprintf(sprintf('\n'))
    toc
