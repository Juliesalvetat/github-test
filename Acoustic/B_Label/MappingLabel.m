% let's map :3
% first things first addpath to m_map
addpath("C:\Users\jsalveta\Desktop\Matlab code\m_map\")
output_path = 'C:\Users\jsalveta\Desktop\Matlab code\figure\';
%load data from EI
filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
Classif = matfile(filename,'Writable',false); clear filename;
Latitude = Classif.Latitude;
Longitude = Classif.Longitude;

% transect map


% Set map bounderies
offset_projection = 0.05; %in degree of Lat and Lon : extend boudaries of map
Longitude_min = -32.62;
Longitude_max = -32.349;
Latitude_min = -3.91;
Latitude_max = -3.79;

figure
%set(gcf, 'Position', get(0, 'Screensize'));
m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
m_line(Longitude,Latitude,'linestyle','n','marker','.','markersize',4,'color',[.5 .5 .5])
m_gshhs_f('patch',[.5 .5 .5]);
m_grid('linestyle','none','box','fancy','fontsize',14,'fontweight','bold');



% load data from Label.mat
path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\Label\Label.mat';
load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;

%Remove Removed Contours
ContourDepth(RemovedContour,:) =[];
ContourIndPing(RemovedContour,:) =[];
LabelIndex(RemovedContour) =[];
ContourIndDepth = NaN(size(ContourDepth));


EsuPingStartEnd = Classif.EsuPingStartEnd;
depth_surface = Classif.depth_surface;
maskFish = Classif.Mask_clean;
boolean_maskFish = ~isnan(maskFish(:,:,1)); %Set boolean matrix echogram ping in x and depth in y axis
Sv_surface = Classif.Sv_surface;
Sv_Fish = Sv_surface.*maskFish;
Sv_Fish = Sv_Fish(:,:,1);
map=colormap(jet(7));

% get CountourDepth index
for i=1:size(ContourDepth,1)
    for j=1:size(ContourDepth,2)
        ind = find(abs(depth_surface-ContourDepth(i,j))==min(abs(depth_surface-ContourDepth(i,j))));
        if ~isempty(ind)
            ContourIndDepth(i,j) = ind(1);
        end
    end
end
count=0;

%construct mask
for i=[1,5,7,8,11,20,21]%3:max(LabelIndex)
    count=count+1;
    BWf = zeros(size(boolean_maskFish));
    
    ContourIndPing_i = ContourIndPing(LabelIndex==i,:);
    ContourIndDepth_i = ContourIndDepth(LabelIndex==i,:);
    
    for j=1:size(ContourIndPing_i,1)
        d=ContourIndPing_i(j,:);
        d=d(~isnan(d));
        for k=1:length(d)
            ind = find(d(k)==EsuPingStartEnd(1,:));
            if ~isempty(ind)
                cip(k)=ind;
            else
                cip(k)=NaN;
            end
        end
        if sum(isnan(cip))<length(cip)-1
            d=cip;d=d(~isnan(cip));
            e=ContourIndDepth_i(j,:);
            e=e(~isnan(e));
            e=e(~isnan(cip));
            
            BW = poly2mask(d,e,size(boolean_maskFish,1),size(boolean_maskFish,2));
            BWf = BWf | BW;
        end
        clear cip
    end
    im_bin = boolean_maskFish & BWf;
    horizontal_sum = sum(im_bin);
    horizontal_sum(horizontal_sum==0)=NaN;
    Lat = Latitude(~isnan(horizontal_sum));
    Lon = Longitude(~isnan(horizontal_sum));
    horizontal_sum = horizontal_sum(~isnan(horizontal_sum));
  if i==1  
    figure
    scatter(Lon,Lat,horizontal_sum,map(count,:)) 
    hold on
  else 
    scatter(Lon,Lat,horizontal_sum,map(count,:))  
  end
end
h1=plot(Lon,Lat,'o')  


Legend=cell(7,1)%  two positions
Legend{1}='bottom.small.fish.school' ;
Legend{2}='small.pelagic.school';
Legend{3}='predators';
Legend{4}='shelfbreak.school';
Legend{5}='shelfbreak.large.fish';
Legend{6}='individual.demersal.fish';
Legend{7}='mixt.reef.fish';
legend(Legend);


    figure
    set(gcf, 'Position', get(0, 'Screensize'));
    m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    m_line(Longitude,Latitude,'linestyle','n','marker','.','markersize',4,'color',[.5 .5 .5])
    m_line(Lon,Lat,'linestyle','n','marker','.','markersize',8,'color','r')
    m_gshhs_f('patch',[.5 .5 .5]);
    m_grid('linestyle','none','box','fancy','fontsize',14,'fontweight','bold');
    title(LabelName{i})
    A = [output_path,LabelName{i},'\FAROFA_2_Label_',LabelName{i},'.png'];
    saveas(gcf,[output_path, 'FAROFA_2_Label_',LabelName{i},'.png'])
    close all
    clear im_bin BWf d e
end


