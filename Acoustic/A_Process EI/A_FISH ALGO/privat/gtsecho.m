function gtsecho(z, y, CLim, Botx)
%gtsecho(z, y, CLim, Botx)
%function to plot an echogram
%IMPUTS:
%   z: Sv matrix
%   y: Depth vector
%   CLIM: Treshold to see Sv values
%   Botx: Depth bottom vector
%OUTPUTS:
%   Echogram
%EXAMPLE
% gtsecho(TData.Sv_conv(:,:,1)',TData.Depth',[-93 -54],TData.Bot)
% load EK500_colourmap.dat;
ek5=ek500();

if nargin==3
    Botx=NaN*ones(length(z));
end
%ek5=ek500();   
%    
%figure   
%set(gcf,'position',get(0,'screensize'))
image(1:size(z,2),y,z,'CDataMapping','scaled');
hold on
plot(Botx,'k')

xlabel('Ping','Fontname','Times New Roman','Fontsize',15)
ylabel('Depth (m)','Fontname','Times New Roman','Fontsize',15)
set(gca,'Fontname','Times New Roman','Fontsize',15)
% xlabel('Ping')
% ylabel('Depth (m)')

colormap(ek5);
caxis(CLim);
colorbar;
set(gca,'Fontname','Times New Roman','Fontsize',15)
hold off


