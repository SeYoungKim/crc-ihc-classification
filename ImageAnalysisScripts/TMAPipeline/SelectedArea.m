%% Display Imgae Functions Adjusted Seg Area to select the regions of interest
%
% Inputs:
%       idx: selected label indexes
%       L: the bw labelled image
%       show: option to show the image. [0,1]= [no, yes]
%       bwlab: option to relabel the new image [0,1]
%       CSeg: the colour original image to draw on
%
%
% keep in L image (0 = false, 1=true)
%
%
% Dependencies: Image Processing Toolbox

function[L, OverOut]=SelectedArea(idx,L, show, bwlab, CSeg)
if ~exist('show','var'), show=0; end
if ~exist('CSeg', 'var'), show=0; end

    se=strel('disk',2);
    Xsize=size(L);

    L2=L;
    a=ismember(reshape(L,Xsize), idx);
    a=imdilate(a, se);
    b=find(a==0);
    L(b)=0;
    
    if show==1 
        L2(b)=255;
        OverOut=SegArea(CSeg, L2); end
    
    if (bwlab==1) L=bwlabel(L,4); end


