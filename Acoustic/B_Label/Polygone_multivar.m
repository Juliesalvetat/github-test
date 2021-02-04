%%%%%$%% Labelling Polygone


addpath('C:\Users\jsalveta\Desktop\Matlab code\gfunctions\')

clear all; clc


for FAROFA=1:3;
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
    count = 0;
    if FAROFA == 1
        path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\Label\Label.mat';
        filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
        
    elseif FAROFA == 2
        path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\Label\Label.mat';
        filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
        
    elseif FAROFA ==3
        path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_3\Label\Label.mat';
        filename = 'C:\Users\jsalveta\Desktop\Data\FAROFA_3\Echointegration_FISH_cleaned\Echointegration.mat';
    end
    % Load Label.mat
    load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;
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
    
    %construct mask
    for i=1:max(LabelIndex)
        
        BWf = zeros(size(boolean_maskFish));
        
        ContourIndPing_i = ContourIndPing(LabelIndex==i,:);
        ContourIndDepth_i = ContourIndDepth(LabelIndex==i,:);
        ContourDepth_i = ContourDepth(LabelIndex==i,:);
        
        for j=1:size(ContourIndPing_i,1)
            disp(['FAROFA ', num2str(FAROFA),' Label : ', num2str(newLabelName{i,2}),' ',num2str(i) ' over ',num2str(max(LabelIndex)), ' labels ; Polygon ',num2str(j),'/', num2str(size(ContourIndPing_i,1))])
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
                Poly_lat = [Poly_lat Latitude(ind_x)];
                Poly_lon = [Poly_lon Longitude(ind_x)];
                Poly_Time = [Poly_Time Time(ind_x)];
                Poly_Bottom = [Poly_Bottom Bottom(ind_x)];
                Poly_depth = [Poly_depth depth_surface(ind_y)];
                Poly_limit_up = [Poly_limit_up min(f)];
                Poly_limit_down = [Poly_limit_down max(f)];
                F = [F FAROFA];
                Poly_number = [Poly_number count];
                Poly_LabelName{count,1} = newLabelName{i,2};
                Poly_EsuPingStart = [Poly_EsuPingStart EsuPingStartEnd(1,ind_x)];
%                 figure; imagesc(BW);xlim([min(d)-100 max(d)+100]);hold on; plot(stat.Centroid(1),stat.Centroid(2),'ro');
                               centroids = cat(1,stat.Centroid);
                               plot(centroids(:,1),centroids(:,2),...
                                   'MarkerEdgeColor','k',...
                                   'MarkerFaceColor',[.49 1 .63],...
                                   'MarkerSize',10)
                %
                BWf = BWf | BW;
                %figure; imagesc(BWf)
                
                %             if length(f)==length(d)
                %                 f=f(~isnan(cip));
                im = boolean_maskFish & BW;
                maskFish = Sv_Fish_70.*im;
                                                 figure;gtsecho(Sv_Fish_70,depth_surface,[-80 -40],Bottom(1,:));
                
                                                 figure;gtsecho(maskFish,depth_surface,[-80 -40],Bottom(1,:)+TransducerDepth);
                                                 hold on
                                                 plot(d,f,'--b','Linewidth',2)
                                                 xlim([min(d)-100 max(d)+100])
                                                 ylim([0 max(max(Bottom(1,min(d)-50:max(d)+50)),max(f))])
                                                 pl=line([ind_x ind_x],[Bottom(ind_x) depth_surface(ind_y)])
                                                 pl.Color = 'black';
                                                 pl.LineStyle = '--';
                                                 pl.LineWidth = 2;
                                                 p2=line([0 ind_x],[depth_surface(ind_y) depth_surface(ind_y)])
                                                 p2.Color = 'black';
                                                 p2.LineStyle = '--';
                                                 p2.LineWidth = 2;                                                 
                                                 plot(ind_x,depth_surface(ind_y),...
                                                     '-mo',...
                                                     'LineWidth',2,...
                                                     'MarkerEdgeColor','k',...
                                                     'MarkerFaceColor',[.49 1 .63],...
                                                     'MarkerSize',10)

                                                
                                                line([0 L],[0 0]);
                %                 saveas(gcf,[path_lib, LabelName2{i},'\FAROFA_2_Label_',num2str(i),'n',num2str(count),'.png'])
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
            Poly_pixel_number = [Poly_pixel_number sum(sum(im_bin))];
            
            Poly_select_Sv_70=NaN(size(im_bin));
            Poly_select_Sa_70=NaN(size(im_bin));
            Poly_select_Sv_200=NaN(size(im_bin));
            Poly_select_Sa_200=NaN(size(im_bin));
            
            Poly_select_Sv_70(im_bin)=Sv_Fish_70(im_bin);
            Poly_select_Sa_70(im_bin)=Sa_Fish_70(im_bin);
            Poly_select_Sv_200(im_bin)=Sv_Fish_200(im_bin);
            Poly_select_Sa_200(im_bin)=Sa_Fish_200(im_bin);
            
            Poly_Sv_70 = Poly_select_Sv_70(~isnan(Poly_select_Sv_70));
            Poly_Sa_70 = Poly_select_Sa_70(~isnan(Poly_select_Sa_70));
            Poly_Sv_200 = Poly_select_Sv_200(~isnan(Poly_select_Sv_200));
            Poly_Sa_200 = Poly_select_Sa_200(~isnan(Poly_select_Sa_200));
            
            Poly_sum_Sa_70 = [Poly_sum_Sa_70 nansum(Poly_Sa_70)];
            Poly_mean_Sv_70 = [Poly_mean_Sv_70 log10(nanmean(10.^(Poly_Sv_70)))];
            Poly_sum_Sa_200 = [Poly_sum_Sa_200 nansum(Poly_Sa_200)];
            Poly_mean_Sv_200 = [Poly_mean_Sv_200 log10(nanmean(10.^(Poly_Sv_200)))];
            
            
            %             figure;gtsecho(D,depth_surface,[-80 -40],Bottom(1,:));
            %             xlim([min(d)-100 max(d)+100])
            
            % figure;imagesc(im_bin);xlim([min(d)-100 max(d)+100])
            maskLabel(im_bin)=i;
            % figure;imagesc(maskLabel);xlim([min(d)-100 max(d)+100]);colorbar
        end
        clear im_bin BWf
    end
    save(['maskLabel_F_',num2str(FAROFA),'.mat'],'maskLabel','Latitude','Longitude','Time','Bottom','-v7.3')
    T = table(Poly_number',Poly_LabelName,Poly_lon',Poly_lat',Poly_Time',Poly_EsuPingStart',Poly_Bottom',...
        Poly_depth',Poly_limit_up',Poly_limit_down',Poly_sum_Sa_70',Poly_sum_Sa_200',Poly_mean_Sv_70',Poly_mean_Sv_200',...
        Poly_pixel_number',F',...
        'VariableNames',{'Poly_LabelName','Poly_number','Poly_Longitude','Poly_Latitude','Poly_Time','Poly_EsuPingStart','Poly_Bottom',...
        'Poly_depth','Poly_limit_up','Poly_limit_down','Poly_sum_Sa_70','Poly_sum_Sa_200','Poly_mean_Sv_70','Poly_mean_Sv_200',...
        'Poly_pixel_number','FAROFA'});
    writetable(T,['F',num2str(FAROFA),'_masklabel.csv'])
    %     figure;imagesc(maskLabel);colorbar
clear    
end
