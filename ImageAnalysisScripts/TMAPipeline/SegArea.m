% Function to highlight selected regions according to a bwlabel object

% Input Required:
%           Samp = original image
%           BWLabO = bw labelled object

function [L2] = SegArea(Samp, BWLabO)
Lrgb = label2rgb(BWLabO, 'jet', 'w', 'shuffle');
Rchan=(0.75*Samp(:,:,1)+0.25*Lrgb(:,:,1));
Bchan=(0.75*Samp(:,:,2)+0.25*Lrgb(:,:,2));
Cchan=(0.75*Samp(:,:,3)+0.25*Lrgb(:,:,3));
L2=cat(3, Rchan, Bchan, Cchan);
