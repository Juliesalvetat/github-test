%%%%%$%% Labelling Polygone
clear all; clc

count = 0;
% Load Label.mat
path2Label = 'C:\Users\jsalveta\Desktop\magiksquare\magiksquare\Cruise_CruiseName\Treatment20200706_134709_FISH\Label\Label.mat';
load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;

% load EI FISH
filename = 'C:\Users\jsalveta\Desktop\magiksquare\DATA\EI FISH\Echointegration.mat';
Classif = matfile(filename,'Writable',false); clear filename;
newLabelName = split(LabelName,'-');

%Remove Removed Contours
ContourDepth(RemovedContour,:) =[];
ContourIndPing(RemovedContour,:) =[];
LabelIndex(RemovedContour) =[];
ContourIndDepth = NaN(size(ContourDepth));

depth_surface = Classif.depth_surface;
% boolean_maskFish = ~isnan(maskFish(:,:,1)); %Set boolean matrix echogram ping in x and depth in y axis
% figure; imagesc(maskFish(:,:,1))
% figure; imagesc(boolean_maskFish)
Time = (Classif.Time / (24*3600))+datenum('01-Jan-1970');
%[Y M D]=datevec(Time);
Bottom = Classif.Bottom(1,:);
TransducerDepth = Classif.TransducerDepth(1,:);
Latitude = Classif.Latitude;
Longitude = Classif.Longitude;
Sv_surface = Classif.Sv_surface;
[~,cent] = min(abs(depth_surface - 100));
maskFish = Classif.Mask_clean(1:cent,:,1);
depth_surface =depth_surface(1:cent);
Sv_Fish_70 = Sv_surface(1:cent,:,1).*maskFish;
Sv_Fish_200 = Sv_surface(1:cent,:,2).*maskFish;
%     figure;gtsecho(Sv_Fish_70,depth_surface,[-80 -40],Bottom(1,:));
%     figure;gtsecho(Sv_Fish(:,:,1),depth_surface,[-80 -40],Bottom(1,:));
%     figure;gtsecho(Sv_surface(:,:,1),depth_surface,[-80 -40],Bottom(1,:));

Sa_surface = Classif.Sa_surface;
Sa_Fish_70 = Sa_surface(1:cent,:,1).*maskFish;
Sa_Fish_200 = Sa_surface(1:cent,:,2).*maskFish;

boolean_maskFish = ~isnan(Sv_Fish_70); %Set boolean matrix echogram ping in x and depth in y axis
EsuPingStartEnd = Classif.EsuPingStartEnd;
clear maskFish %Sv_surface
% get CountourDepth index
for i=1:size(ContourDepth,1)
    for j=1:size(ContourDepth,2)
        ind = find(abs(depth_surface-ContourDepth(i,j))==min(abs(depth_surface-ContourDepth(i,j))));
        if ~isempty(ind)
            ContourIndDepth(i,j) = ind(1);
        end
    end
end
maskLabel = zeros(size(boolean_maskFish));
%path_lib = 'C:\Users\jsalveta\Desktop\Matlab code\Labelling Code\savingMaskforAMD Library\';
Poly_lat = [];
    Poly_lon = [];
    Poly_Time = [];
    Poly_Bottom = [];
    Poly_depth = [];
    Poly_sum_Sa_70 = [];
    Poly_sum_Sa_200 = [];
    Poly_mean_Sv_70 = [];
    Poly_mean_Sv_200 = [];
    Poly_pixel_number = [];
    F = [];
    Poly_number = [];
    Poly_LabelName =  strings(1);
    Poly_limit_up = [];
    Poly_limit_down = [];
    Poly_EsuPingStart = [];
%construct mask
for i=1:max(LabelIndex)
    
    BWf = zeros(size(boolean_maskFish));
    
    ContourIndPing_i = ContourIndPing(LabelIndex==i,:);
    ContourIndDepth_i = ContourIndDepth(LabelIndex==i,:);
    ContourDepth_i = ContourDepth(LabelIndex==i,:);
    
    for j=1:size(ContourIndPing_i,1)
        disp([' Label : ', num2str(newLabelName{i,2}),' ',num2str(i) ' over ',num2str(max(LabelIndex)), ' labels ; Polygon ',num2str(j),'/', num2str(size(ContourIndPing_i,1))])
        %t=sprintf('lig1\nlig2');
        c=ContourIndPing_i(j,:);
        c=c(~isnan(c));
        cip = NaN(1,length(c));
        for k=1:length(c)
            [~,ind]=min(abs(c(k)-EsuPingStartEnd(1,:)));
            %ind = find(c(k)==EsuPingStartEnd(1,:));
            if ~isempty(ind)
                cip(k)=ind;
            end
        end
        
        if sum(isnan(cip))<length(cip)-1
            count=count+1;
            d=cip;d=d(~isnan(cip));
            e=ContourIndDepth_i(j,:);
            e=e(~isnan(e));
            e=e(~isnan(cip));
            f=ContourDepth_i(j,:);
            f=f(~isnan(f));
            BW = poly2mask(d,e,size(boolean_maskFish,1),size(boolean_maskFish,2));
            stat = regionprops(BW ,'Centroid','Area');
            area = cat(1,stat.Area);
            if length(area)>1
                centroids = cat(1,stat.Centroid);
                [a m] = max(area);
                ind_x = round(centroids(m,1));
                ind_y = round(centroids(m,2));
            else
                ind_x = round(stat.Centroid(1));
                ind_y = round(stat.Centroid(2));
            end
            Poly_lat(count) = Latitude(ind_x);
            Poly_lon(count) =  Longitude(ind_x);
            Poly_Time(count) = Time(ind_x);
            Poly_Bottom(count) = Bottom(ind_x);
            Poly_depth(count) = depth_surface(ind_y);
            Poly_limit_up(count) = min(f);
            Poly_limit_down(count) = max(f);
            Poly_number(count) = count;
            Poly_LabelName{count,1} = newLabelName{i,2};
            Poly_EsuPingStart(count) =  EsuPingStartEnd(1,ind_x);
            %                 figure; imagesc(BW);xlim([min(d)-100 max(d)+100]);hold on; plot(stat.Centroid(1),stat.Centroid(2),'ro');
            %                                centroids = cat(1,stat.Centroid);
            %                 plot(centroids(:,1),centroids(:,2),'b*')
            %
            BWf = BWf | BW;
            %figure; imagesc(BWf)
            
            %             if length(f)==length(d)
            %                 f=f(~isnan(cip));
            im = boolean_maskFish & BW;
            maskFish = Sv_Fish_70.*im;
            %                                 figure;gtsecho(Sv_Fish_70,depth_surface,[-80 -40],Bottom(1,:));
            
            %                                 figure;gtsecho(maskFish,depth_surface,[-80 -40],Bottom(1,:)+TransducerDepth);
            %                                 hold on
            %                                 plot(d,f,'--b','Linewidth',2)
            %                                 xlim([min(d)-100 max(d)+100])
            %                                 ylim([0 max(max(Bottom(1,min(d)-50:max(d)+50)),max(f))])
            % %                 saveas(gcf,[path_lib, LabelName2{i},'\FAROFA_2_Label_',num2str(i),'n',num2str(count),'.png'])
            % %                 close all
            %                 % figure;imagesc(BW);xlim([min(d)-100 max(d)+100])
            %                 % figure;imagesc(boolean_maskFish);xlim([min(d)-100 max(d)+100])
            %                 % figure;imagesc(Sv_Fish);xlim([min(d)-100 max(d)+100])
            %                 clear d e f cip im maskFish
            % %             end
            %         else
            % figure;imagesc(BWf);xlim([min(d)-100 max(d)+100]);colorbar
            %                 figure;imagesc(im_bin);xlim([min(d)-100 max(d)+100]);colorbar
            %                figure;imagesc(Sa_Fish_70(im_bin));xlim([min(d)-100 max(d)+100]);colorbar
            %                figure;imagesc(Sa_Fish_70);xlim([min(d)-100 max(d)+100]);colorbar; caxis([0 40])
            %                figure;imagesc(Sa_Fish_200);xlim([min(d)-100 max(d)+100]);colorbar; caxis([0 40])
            %                figure;imagesc(Poly_select_Sa_200);xlim([min(d)-100 max(d)+100]);colorbar; caxis([0 40])
            %                figure;imagesc(Poly_select_Sv_200);xlim([min(d)-100 max(d)+100]);colorbar;
            %                figure;imagesc(Poly_select_Sv_70);xlim([min(d)-100 max(d)+100]);colorbar;
            %                figure;gtsecho(Poly_select_Sa_70,depth_surface,[0 40],Bottom(1,:));xlim([min(d)-100 max(d)+100])
            
            clear d e f cip
        end
        im_bin = boolean_maskFish & BWf;
        Poly_pixel_number(count) = sum(sum(im_bin));
        
        Poly_select_Sv_70=NaN(size(im_bin));
        Poly_select_Sa_70=NaN(size(im_bin));
        Poly_select_Sv_200=NaN(size(im_bin));
        Poly_select_Sa_200=NaN(size(im_bin));
        
        Poly_select_Sv_70(im_bin)=Sv_Fish_70(im_bin);
        Poly_select_Sa_70(im_bin)=Sa_Fish_70(im_bin);
        Poly_select_Sv_200(im_bin)=Sv_Fish_200(im_bin);
        Poly_select_Sa_200(im_bin)=Sa_Fish_200(im_bin);
        
            Poly_Sv_70 = Poly_select_Sv_70(~isnan(Poly_select_Sv_70));
            Poly_Sv_200 = Poly_select_Sv_200(~isnan(Poly_select_Sv_200));
            tmp = sum(~isnan(Poly_select_Sa_70),1);
            nb_ping_non_nul = length(find(tmp>0));
            
            Poly_Sa_70 = Poly_select_Sa_70(~isnan(Poly_select_Sa_70));
            Poly_Sa_200 = Poly_select_Sa_200(~isnan(Poly_select_Sa_200));
            
            Poly_sum_Sa_70 = [Poly_sum_Sa_70 nansum(Poly_Sa_70)/nb_ping_non_nul];
            Poly_mean_Sv_70 = [Poly_mean_Sv_70 log10(nanmean(10.^(Poly_Sv_70)))];
            Poly_sum_Sa_200 = [Poly_sum_Sa_200 nansum(Poly_Sa_200)/nb_ping_non_nul];
            Poly_mean_Sv_200 = [Poly_mean_Sv_200 log10(nanmean(10.^(Poly_Sv_200)))];
            
        
        
        %             figure;gtsecho(D,depth_surface,[-80 -40],Bottom(1,:));
        %             xlim([min(d)-100 max(d)+100])
        
        % figure;imagesc(im_bin);xlim([min(d)-100 max(d)+100])
        maskLabel(im_bin)=i;
        % figure;imagesc(maskLabel);xlim([min(d)-100 max(d)+100]);colorbar
    end
    clear im_bin BWf
end
save('PolyLabel_Cris.mat','maskLabel','Latitude','Longitude','Time','Bottom','-v7.3')
T = table(Poly_number',Poly_LabelName,Poly_lon',Poly_lat',Poly_Time',Poly_EsuPingStart',Poly_Bottom',...
    Poly_depth',Poly_limit_up',Poly_limit_down',Poly_sum_Sa_70',Poly_sum_Sa_200',Poly_mean_Sv_70',Poly_mean_Sv_200',...
    Poly_pixel_number',...
    'VariableNames',{'Poly_LabelName','Poly_number','Poly_Longitude','Poly_Latitude','Poly_Time','Poly_EsuPingStart','Poly_Bottom',...
    'Poly_depth','Poly_limit_up','Poly_limit_down','Poly_sum_Sa_70','Poly_sum_Sa_200','Poly_mean_Sv_70','Poly_mean_Sv_200',...
    'Poly_pixel_number'});
writetable(T,['PolyLabel_Cris.csv'])
%     figure;imagesc(maskLabel);colorbar
clear

clear all;close all;clc

%%%%% Add some stuff to Polygone_multivar table
% shelfbreak distance
% coast distance

% add m_map path
addpath('C:\Users\jsalveta\Desktop\Matlab code\m_map\')
% set boundarie limits
Longitude_min = -32.5; % (without Drina submount)
Longitude_max = -32.349;
Latitude_min = -3.91;
Latitude_max = -3.79;
offset_projection = 0.005;


load('C:\Users\jsalveta\Desktop\Matlab code\bathy_coast_FdN\poubelle\Bati_Ilhas.mat')
offset_projection = 0.05;
ind_lon = find(lon > Longitude_min - offset_projection & lon < Longitude_max + offset_projection);
ind_lat = find(lat > Latitude_min - offset_projection & lat < Latitude_max + offset_projection);
Elevation =Elevation(ind_lat,ind_lon);
%figure
[C,h] = contour(lon(ind_lon),lat(ind_lat), Elevation,[-80 -80]); %hold on
Lvls = h.LevelList;
idxc = find(C(1,:) == Lvls);
Llen = C(2,idxc);
conturc{1} = C(:,idxc(1)+1 : idxc(1)+1+Llen(1)-1);
conturc{2} = C(:,idxc(2)+1 : idxc(2)+1+Llen(2)-1);
% figure
% plot(conturc{1}(1,:),conturc{1}(2,:))
% hold on
% plot(conturc{2}(1,:),conturc{2}(2,:),'r')
ShelfbreakLong = conturc{1}(1,:);
ShelfbreakLat = conturc{1}(2,:);


% read Polygone table from Polygone_multivar.m
Poly_table = readtable('PolyLabel_Cris.csv');
Longitude = Poly_table.Poly_Longitude;
Latitude = Poly_table.Poly_Latitude;
Bottom = Poly_table.Poly_Bottom;


%%%%%%%%%%% distance to the coast
%%% small twist coast is actually FN
load('C:\Users\jsalveta\Desktop\Matlab code\bathy_coast_FdN\FN.mat')
long = ncst(:,1);
lat = ncst(:,2);

coastLong = long(~isnan(long));
coastLat = lat(~isnan(lat));

lat = Latitude;        %# Inland latitude points (in degrees)
long = Longitude;  %# Inland longitude points (in degrees)
nPoints = numel(lat);         %# Number of map points
scale = pi/180;               %# Scale to convert degrees to radians
radiusEarth = 3958.76;        %# Average radius of Earth, in miles
distanceToCoast = zeros(1,nPoints);   %# Preallocate distance measure
coastIndex = zeros(1,nPoints);        %# Preallocate a coastal point index
for iPoint = 1:nPoints                %# Loop over map points
    rho = cos(scale.*lat(iPoint)).*...  %# Compute central angles from map
        cos(scale.*coastLat).*...     %#   point to all coastal points
        cos(scale.*(coastLong-long(iPoint)))+...
        sin(scale.*lat(iPoint)).*...
        sin(scale.*coastLat);
    d = radiusEarth.*acos(rho);         %# Compute great-circle distances
    [distanceToCoast(iPoint),coastIndex(iPoint)] = min(d);  %# Find minimum
end
distanceToCoast = distanceToCoast'.*1.60934;
% visual check
figure
plot(coastLong,coastLat,'.');hold on
for i=1:length(long)
plot([long(i)  coastLong(coastIndex(i))],...    %'# Plot the inland points and
     [lat(i) coastLat(coastIndex(i))],...      %'#   nearest coastal points
     'ro-');
end




% str = strcat(num2str(distanceToCoast(i),...  %'# Make text for the distances
%                      '%0.1f'),{' km'});
% text(long(i),lat(i),str(i),'Color','r','VerticalAlignment','bottom');  %# Plot the text



%%%%%%%%%%%% distance au talu 80m
load('C:\Users\jsalveta\Desktop\Matlab code\bathy_coast_FdN\poubelle\Bati_Ilhas.mat')
Longitude_min = -32.5; % (without Drina submount)
Longitude_max = -32.349;
Latitude_min = -3.91;
Latitude_max = -3.79;
offset_projection = 0.05;
ind_lon = find(lon > Longitude_min - offset_projection & lon < Longitude_max + offset_projection);
ind_lat = find(lat > Latitude_min - offset_projection & lat < Latitude_max + offset_projection);
Elevation =Elevation(ind_lat,ind_lon);


figure
[C,h] = contour(lon(ind_lon),lat(ind_lat), Elevation,[-80 -80]); %hold on
Lvls = h.LevelList;
    idxc = find(C(1,:) == Lvls);
    Llen = C(2,idxc);
    conturc{1} = C(:,idxc(1)+1 : idxc(1)+1+Llen(1)-1);
    conturc{2} = C(:,idxc(2)+1 : idxc(2)+1+Llen(2)-1);
% figure
% plot(conturc{1}(1,:),conturc{1}(2,:))
% hold on
% plot(conturc{2}(1,:),conturc{2}(2,:),'r')
%
%
% ShelfbreakLong = conturc{1}(1,:);
% ShelfbreakLat = conturc{1}(2,:);

lat = Latitude;        %# Inland latitude points (in degrees)
long = Longitude;  %# Inland longitude points (in degrees)
nPoints = numel(lat);         %# Number of map points
scale = pi/180;               %# Scale to convert degrees to radians
radiusEarth = 3958.76;        %# Average radius of Earth, in miles
distanceToShelfbreak = zeros(1,nPoints);   %# Preallocate distance measure
ShelfbreakIndex = zeros(1,nPoints);        %# Preallocate a coastal point index
for iPoint = 1:nPoints                %# Loop over map points
    rho = cos(scale.*lat(iPoint)).*...  %# Compute central angles from map
        cos(scale.*ShelfbreakLat).*...     %#   point to all coastal points
        cos(scale.*(ShelfbreakLong-long(iPoint)))+...
        sin(scale.*lat(iPoint)).*...
        sin(scale.*ShelfbreakLat);
    d = radiusEarth.*acos(rho);         %# Compute great-circle distances
    [distanceToShelfbreak(iPoint),ShelfbreakIndex(iPoint)] = min(d);  %# Find minimum
end
distanceToShelfbreak = distanceToShelfbreak'.*1.60934;
% visual check
    figure
plot(coastLong,coastLat,'.');hold on
    plot(ShelfbreakLong,ShelfbreakLat,'-');
    for i=1:length(long)
        plot([long(i)  ShelfbreakLong(ShelfbreakIndex(i))],...    %'# Plot the inland points and
            [lat(i) ShelfbreakLat(ShelfbreakIndex(i))],...      %'#   nearest coastal points
            'ro-');
    end

Poly_table.distanceToCoast = distanceToCoast;
Poly_table.distanceToShelfbreak = distanceToShelfbreak;
writetable(Poly_table,'PolyLabel_Cris.csv')

