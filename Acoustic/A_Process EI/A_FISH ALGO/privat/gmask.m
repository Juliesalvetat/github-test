function datai=gmask(data,mask,oval,inv)
%datai=gmask(data,mask,oval,inv)
%function to apply a boolean mask
%IMPUTS:
%   data: Sv matrix
%   mask: logical matrix
%   oval: Value to be changed to false case
%   inv: just if you want to invert the logical matrix
%OUTPUTS:
%   Sv matrix
%
%EXAMPLE
% S1.mask=Sv_conv(:,:,1)>-60;
% S1.Sv=gmask(Sv_conv(:,:,1),S1.mask,-999);

if nargin==3
    inv=0;
end
datai=oval*ones(size(data));
if inv==1
    mask=~mask;
end
    datai(mask)=data(mask);
    datai(isnan(data))=NaN;
