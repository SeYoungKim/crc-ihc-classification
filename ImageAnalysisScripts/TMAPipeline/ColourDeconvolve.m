function [imOut]=ColourDeconvolve(imInO, Outtype, D)
% COLOURDECONVOLVE: Deconvolves TMA image into DAB, E & H components 
%   
%   imOut = ColourDeconvolve(imIn) creates an output image (imOut) containing a
%   projection of the RGB input colour image (imIn) in three dimensions:
%   (Haematoxylin, Eosin and DAB), using the deconvolution matrix specified
%   by Ruifrok(2001).
%   I.e. 
%       imOut(:,:,1) is the Haematoxylin channel
%       imOut(:,:,2) is the Eosin channel
%       imOut(:,:,3) is the DAB (Brown) channel
%
%   imOut = ColourDeconvolve(imIn, D) deconvolves the image according to a
%   (3 colour x 3 stain) matrix D, containing experimental absorbances in 
%   the RGB channels for 3 stains of interest (H, E, DAB). 
%       Outtype - indicates the colour conversion. Either 'OD' or 'RGB'
%                 Default to 'OD' ie. from RGB -->OD
%
%   Example:
%   --------
%       TMAOD = ColourDeconvolve(TMA1, 'OD');
%       TMAlight = ColourDeconvolve(TMAOD, 'RGB');
%       imtool(TMAOD);
%       imtool(TMAlight);
%
%   Reference:
%   ---------
%   Ruifrok AC (2001), Quantification of histochemical staining by color 
%   deconvolution, Analytical and Quantitative Cytology and Staining
%   23(4):291-9

%   Check inputs
tic

    N = nargin;
    error(nargchk(1,2,N));
    assert(ndims(imInO)==3, 'Dimensions of image incorrect. Enter a RGB image.')
    if exist('Outtype', 'var')
        assert(strcmp(Outtype, 'OD')==1|strcmp(Outtype, 'RGB')==1, 'Indicate either ''OD'' or ''RGB'' for image transform');
    else
        Outtype='OD';
    end
        nn=size(imInO);
    
    
%   Specify Ruifrok's deconvolution matrix if there is no input
    if ~exist('D', 'var')
        switch Outtype
            case 'OD'
            Dinv=[ 1.9382   -1.0449   -0.6056;   -0.1034    1.1210   -0.0721; -0.5934   -0.4514    1.5430]';
            case 'RGB'
            Dinv=[0.6412    0.7124    0.2849;   0.0764    0.9941    0.0765; 0.2689    0.5648    0.7800]'; 
        end
     else 
        assert(rank(D)==3, 'Input matrix D not linearly indepedent!')
        Dinv=D;
    end

% Convert to OD space if specified    
    switch Outtype
        case 'OD'
             imIn=-log10(im2double(imInO-1)+0.01);
        case 'RGB'
             imIn=im2double(imInO);
    end
  
%   separate colours into RGB and rotate matrix
    R=LineariseChannel(imIn(:,:,1));
    G=LineariseChannel(imIn(:,:,2));
    B=LineariseChannel(imIn(:,:,3));
    NewMat=[R';G';B'];
%   Deconvolute
    DeconMatrix=Dinv*NewMat;
%   Convert the Arrays into Matrices
    H=DeconMatrix(1,:);
    E=DeconMatrix(2,:);
    DAB=DeconMatrix(3,:);
    H=reshape(H, nn(1:2));
    E=reshape(E, nn(1:2));
    DAB=reshape(DAB, nn(1:2));
%   Colour white space in original image
    x1=find(imInO(:,:,1)==255 & imInO(:,:,2)==255 & imInO(:,:,3)==255);    
    E(x1)=0;
    DAB(x1)=0;
    H(x1)=0;
%   Concatenate to obtain new image
    imOut=cat(3,H,E,DAB);
    if strcmp(Outtype,'RGB')==1
        imOut=10.^(-imOut);
    end
    imOut=im2uint8(imOut);
     fprintf(sprintf('\n'))
    toc
end


function [Channel]=LineariseChannel(Channel)
% Linearise a channel and convert from uint8 to double
    Channel=Channel(:);
  end

% function [Channel]=MatrixChannel(Channel, size)
% % Convert channel from double to uint8 and reshape matrix
%      Channel=reshape(Channel, size);
% end