%%%% Average Horizontal Sa & calculate rugosity
% Use Echointegration.mat and Sa_Horizontal.mat produced by
% A_Calculate_Sa_horizontal.mat
% Does it on every folder of each farofa cruise

clc; clear all; close all


% Select and load Echointegration file (FISH or NO FISH)
% to load CumulatedGPSDistanceMeter, Time, DateTime, Longitude,
% Latitude, Bottom, EsuPingStartEnd
EI_path = ['.\FISH\Echointegration.mat'];
EI = matfile(EI_path,'Writable',false);

% load Horizontal_sa_surface 1 ping fish and no fish
SaHorizontal_path_fish = ['.\FISH\Sa_Horizontal.mat'];
SaHorizontal_fish = matfile(SaHorizontal_path_fish,'Writable',false);
SaHorizontal_path_nofish = ['.\NO FISH\Sa_Horizontal.mat'];
SaHorizontal_nofish = matfile(SaHorizontal_path_nofish,'Writable',false);

% Interval to apply in meters
Interval=25;

%load variable of interest
CumulatedGPSDistanceMeter = EI.CumulatedGPSDistanceMeter;
DateTime = EI.DateTime;
Longitude = EI.Longitude;
Latitude = EI.Latitude;
Time = (EI.Time / (24*3600))+datenum('01-Jan-1970');
Bottom = EI.Bottom;
EsuPingStartEnd = EI.EsuPingStartEnd;
Horizontal_sa_surface_fish = SaHorizontal_fish.Horizontal_sa_surface;
Horizontal_sa_surface_nofish = SaHorizontal_nofish.Horizontal_sa_surface;

% check Sa fish /no fish
%     figure
%     subplot(3,1,[1:2]);plot(Horizontal_sa_surface_fish(:,1));hold on;plot(Horizontal_sa_surface_nofish(:,1))
%     title('SaSum Fish (blue) no fish (red)')
%     ylabel('Sa 1 ping')
%     subplot(3,1,3);plot(-Bottom(1,:))
%     ylabel('Bottom')
%     xlabel('Pings')
%

% defining interval indexes
dist_index=ceil(CumulatedGPSDistanceMeter/Interval)*Interval;
fix_index=unique(dist_index);
dist_mat=min(fix_index):Interval:max(fix_index);

% Calculating
for i=2:length(fix_index)
    int_ind_EI=find(dist_index==fix_index(i));
    int_pos_table=find(dist_mat==fix_index(i));
 


%if  min(dist_mat(interval_position) - CumulatedGPSDistanceMeter(interval_indexes))<Interval*0.1
%figure;plot(fix_index)
interval_vector_char(int_pos_table) = string([num2str(fix_index(i-1)),'-',num2str(fix_index(i))]);
DATE_start(int_pos_table) = string(DateTime(int_ind_EI(1),1:10));
DATE_end(int_pos_table) = string(DateTime(int_ind_EI(end),1:10));
interval_edge_diff(int_pos_table) = min(fix_index(i) - CumulatedGPSDistanceMeter(int_ind_EI));
ESU_start(int_pos_table) = EsuPingStartEnd(1,int_ind_EI(1));
ESU_end(int_pos_table) = EsuPingStartEnd(2,int_ind_EI(end));
GPSDistanceMeter(int_pos_table) = CumulatedGPSDistanceMeter(int_ind_EI(end));
nb_pings_by_interval(int_pos_table) = length(int_ind_EI);
Time_start(int_pos_table) = Time(int_ind_EI(1));
Time_end(int_pos_table) = Time(int_ind_EI(end));
Lat_start(int_pos_table) = Latitude(int_ind_EI(1));
Lat_end(int_pos_table) = Latitude(int_ind_EI(end));
Lon_start(int_pos_table) = Longitude(int_ind_EI(1));
Lon_end(int_pos_table) = Longitude(int_ind_EI(end));
Depth_start(int_pos_table) = Bottom(1,int_ind_EI(1));
Depth_end(int_pos_table) = Bottom(1,int_ind_EI(end));
Depth_std(int_pos_table) = nanstd(Bottom(1,int_ind_EI));
Time_mean(int_pos_table) = nanmean(Time(int_ind_EI));
Lat_mean(int_pos_table) = nanmean(Latitude(int_ind_EI));
Lon_mean(int_pos_table) = nanmean(Longitude(int_ind_EI));
Depth_mean(1,int_pos_table)= nanmean(Bottom(1,int_ind_EI));
SaFish(int_pos_table,:) = nanmean(Horizontal_sa_surface_fish(int_ind_EI,:),1);
SaNoFish(int_pos_table,:) = nanmean(Horizontal_sa_surface_nofish(int_ind_EI,:),1);

% rugositee
% distance entre les points
temp_cumul = CumulatedGPSDistanceMeter(int_ind_EI);
temp_depth = Bottom(1,int_ind_EI);
slope(int_pos_table) = (temp_depth(end)-temp_depth(1))/(temp_cumul(end)-temp_cumul(1)); %pente

for r = 1:length(temp_cumul)-1
    distance(r) = sqrt((temp_cumul(r+1)-temp_cumul(r))^2 + (temp_depth(r+1)-temp_depth(r))^2); %distance
end
rug(int_pos_table) = 1 - ((temp_cumul(end)-temp_cumul(1))/sum(distance));

end


%check rugosity
%     figure
%     subplot(3,1,1);plot(rug)
%     ylabel('Rugosity')
%     subplot(3,1,[2:3]);plot(-Depth_mean)
%     ylabel('Bottom')
%     xlabel('Pings')
%     end


%     for i=1:length(SaNoFish)
%         DATE(i+ day_start_ind(d)-1) = string(DateTime(day_start_ind(d),1:10));
%         interval_vector_char(i+ day_start_ind(d)-1) = string([num2str(dist_mat(i-1)),'-',num2str(dist_mat(i))]);
%         F(i+ day_start_ind(d)-1) = FAROFA;
%
%     end
%figure; plot(interval_edge_diff)
% Saving data
save('Averaged_SaSum_magic_square.mat','DATE_start','DATE_end','interval_vector_char','GPSDistanceMeter','interval_edge_diff','ESU_start','ESU_end',...
    'nb_pings_by_interval','Time_start','Time_end','Lat_start','Lat_end','Lon_start','Lon_end',...
    'Depth_start','Depth_end','Depth_std','rug','slope','Time_mean','Lat_mean','Lon_mean','Depth_mean','SaFish','SaNoFish','-v7.3')

SaFish_70 = SaFish(:,1);
SaFish_200 = SaFish(:,2);
SaNoFish_70 = SaNoFish(:,1);
SaNoFish_200 = SaNoFish(:,2);


T = table(DATE_start',DATE_end',interval_vector_char',interval_edge_diff',GPSDistanceMeter',ESU_start',ESU_end',nb_pings_by_interval',Time_start',Time_end',...
    Lat_start',Lat_end',Lon_start',Lon_end',Depth_start',Depth_end',Depth_std',rug',slope',Time_mean',Lat_mean',Lon_mean',...
    Depth_mean',SaFish_70,SaFish_200,SaNoFish_70,SaNoFish_200,...
    'VariableNames',{'DATE_start','DATE_end','interval_vector_char','interval_edge_diff','GPSDistanceMeter','ESU_start','ESU_end',...
    'nb_pings_by_interval','Time_start','Time_end','Lat_start','Lat_end','Lon_start','Lon_end',...
    'Depth_start','Depth_end','Depth_std','Rugosity','Slope','Time_mean','Lat_mean','Lon_mean','Depth_mean','SaFish_70','SaFish_200',...
    'SaNoFish_70','SaNoFish_200'});

% Empty and remove rows : non existing interval & high mismatch with
% interval edges
% matlab automatically fill non existing space (non existing intervals) in vector with zeros
ind_empty_interval = find(nb_pings_by_interval==0);
repNaN = cell(length(ind_empty_interval), width(T));
repNaN(:) = {NaN};
T(ind_empty_interval,:) = repNaN;

%big differences betwen edges of interval and real data
ind_big_edge_diff = find(interval_edge_diff>Interval*0.2); % 20% 5m for 25m interval
for i=1:length(ind_big_edge_diff)
    ind = find(ESU_start(ind_big_edge_diff(i)) == EsuPingStartEnd(1,:));
    IND(i) =ind(1);
end
% %     figure;plot(interval_edge_diff(ind_big_edge_diff))
% %         figure;plot(GPSDistanceMeter);hold on;plot(GPSDistanceMeter(ind_big_edge_diff),'*r')
% %         plot(repmat(xlim',1,size(IND,2)),[GPSDistanceMeter(ind_big_edge_diff); GPSDistanceMeter(ind_big_edge_diff)], '-r')
% %     
%     figure;
%     plot(CumulatedGPSDistanceMeter,'.'); hold on
%     %plot(repmat(xlim',1,size(IND,2)),[GPSDistanceMeter(ind_big_edge_diff); GPSDistanceMeter(ind_big_edge_diff)], '-r')
%     plot(IND,GPSDistanceMeter(ind_big_edge_diff),'*r')
%     xlabel('Pings')
%     ylabel('CumulatedGPSDistanceMeter')
%     title(['Interval edge mismatch >',num2str(Interval*0.2),' m ; ',num2str(length(ind_big_edge_diff)),' red stars '])
%
%     plot(repmat(xlim',1,size(compare_vec,2)),[compare_vec; compare_vec], '-r')
%     ylabel('CumulatedGPSDistanceMeter')
%     title(['day ',num2str(i),' : ', DateTime(I(i),1:10)])

repNaN = cell(length(ind_big_edge_diff), width(T));
repNaN(:) = {NaN};
T(ind_big_edge_diff,:) = repNaN;
% remove NaN rows
T = rmmissing(T);

% save in csv
writetable(T,'magic_square_Averaged_SaSum.csv')
clear all

