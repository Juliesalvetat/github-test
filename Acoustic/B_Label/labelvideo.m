%%%%%
%is this label inside a video?
%this will save all label that are part of a video or 5min after and before
%mainly for blob because I forgot to implement a blob.vid label or
%bottom.small.fish.school.vid
addpath(".\gfunctions")
    path_lib = 'C:\Users\jsalveta\Desktop\Matlab code\Labelling\blob.vid\';

% You need video Times for FAROFA2 from logbook
load('C:\Users\jsalveta\Desktop\processing video\data for AddTime\DateTime_VP_DV_FAROFA_2.mat')
%you can joinVP and DVmaybe that interresting? be carrefull might contain
%datetime from empty videos
DateTimeStart=[DateTimeStartVP;DateTimeStartDV];
DateTimeEnd=[DateTimeEndVP;DateTimeEndDV];
% ransform into dateTime format
DateTimeStart = datetime(DateTimeStart,'InputFormat','dd/MM/yyyy HH:mm');
DateTimeEnd = datetime(DateTimeEnd,'InputFormat','dd/MM/yyyy HH:mm');

% and the Time in your EI
filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
Classif = matfile(filename,'Writable',false); clear filename;
DateTime = Classif.DateTime;
%change DateTime format char to datetime
t = datetime(DateTime(:,1:16),'InputFormat','yyyy-MM-dd HH:mm');

% you also need Sv matrice for ploting
depth_surface = Classif.depth_surface;
Bottom = Classif.Bottom;
Sv_Fish = Classif.Sv_surface(:,:,1);
EsuPingStartEnd = Classif.EsuPingStartEnd;

%funtion isbetween to determinate if the label is in between the date and time of
%video
% here we want to plot each echogram corresponding to a video ideally I
% will increament 5min before and 5 minute after with a lkne to mark the
% start and end f the video
for i=1:size(DateTimeStart,1)
    tlower=DateTimeStart(i,:);
    tupper=DateTimeEnd(i,:);
    tf = isbetween(t,tlower,tupper);
    ind_section = find(tf==1);
    section = min(ind_section):max(ind_section);
    Sv_Fish_section = Sv_Fish(:,section);
    gtsecho(Sv_Fish_section,depth_surface,[-80 -40],Bottom(1,section));
    
    if i<35
        saveas(gcf,['C:\Users\jsalveta\Desktop\Matlab code\Labelling\video\FAROFA_2_VP_',num2str(i),'.png'])
    else
        saveas(gcf,['C:\Users\jsalveta\Desktop\Matlab code\Labelling\video\FAROFA_2_DV_',num2str(i-34),'.png'])
        %saveas(gcf,[path_lib,'\FAROFA_2_video_n',num2str(count),'.png'])
    end
    
    %don't know yet what name to use
    %saveas(gcf,[path_lib,'\FAROFA_2_video_n',num2str(count),'.png'])
    
end
% Now we include label countours

% Load Label.mat
path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\Label\Label.mat';
load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;

%Remove Removed Contours
ContourDepth(RemovedContour,:) =[];
ContourIndPing(RemovedContour,:) =[];
LabelIndex(RemovedContour) =[];
ContourIndDepth = NaN(size(ContourDepth));
% ContourIndPingEI = NaN(size(ContourDepth));

% get CountourDepth index
for i=1:size(ContourDepth,1)
    for j=1:size(ContourDepth,2)
        ind = find(abs(depth_surface-ContourDepth(i,j))==min(abs(depth_surface-ContourDepth(i,j))));
        if ~isempty(ind)
            ContourIndDepth(i,j) = ind(1);
        end
    end
end

% firts lets try only for bottom.small.fish.school i=1
ContourIndPingESU = NaN(size( ContourIndPing));

for i=1:size(ContourIndPing,1)
    for j=1:size(ContourIndPing,2)
        ind = find (ContourIndPing(i,j)==EsuPingStartEnd(1,:));
        if ~isempty(ind)
            ContourIndPingESU(i,j)=ind;
        else
            ContourIndPingESU(i,j)=NaN;
            
        end
    end
end

for i=1:max(LabelIndex)
    
    ContourIndPing_i = ContourIndPingESU(LabelIndex==i,:);
    ContourIndDepth_i = ContourIndDepth(LabelIndex==i,:);
    ContourDepth_i = ContourDepth(LabelIndex==i,:);
    
    for j=1:size(ContourIndPing_i,1)
        c=ContourIndPing_i(j,:);
        d=c(~isnan(c));
        
        for k=1:size(DateTimeStart,1)
            tlower=DateTimeStart(k,:);
            tupper=DateTimeEnd(k,:);
            tf = isbetween(t(d),tlower,tupper);
            if sum(tf)>=1                
                e=ContourIndDepth_i(j,:);
                e=e(~isnan(c));
                f=ContourDepth_i(j,:);
                f=f(~isnan(c));
                figure
                gtsecho(Sv_Fish,depth_surface,[-80 -40],Bottom(1,:));
                hold on
                plot(d,f,'--b','Linewidth',2)
                xlim([min(d)-100 max(d)+100])
                ylim([0 max(max(Bottom(1,min(d)-50:max(d)+50)),max(f))])
                if k<35
                saveas(gcf,[path_lib, 'FAROFA_2_VP_',num2str(k),'.png'])
                else                   
                saveas(gcf,[path_lib, 'FAROFA_2_DV_',num2str(k-34),'.png'])
                end
               
            end
        end
    end
end
min(ContourIndPing_i)
ind = find(d(k)==EsuPingStartEnd(1,:));

for j=1:size(ContourIndPing_i,1)
    count = count +1;
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
        f=ContourDepth_i(j,:);
        f=f(~isnan(f));
        f=f(~isnan(cip));
        BW = poly2mask(d,e,size(boolean_maskFish,1),size(boolean_maskFish,2));
        im_bin = boolean_maskFish & BW;
        maskLabel = Sv_Fish.*im_bin;
        gtsecho(maskLabel,depth_surface,[-80 -40],Bottom(1,:));
        hold on
        plot(d,f,'--b','Linewidth',2)
        xlim([min(d)-100 max(d)+100])
        ylim([0 max(max(Bottom(1,min(d)-50:max(d)+50)),max(f))])
        saveas(gcf,[path_lib, LabelName{i},'\FAROFA_2_Label_n',num2str(count),'.png'])
        close all
        clear cip
    else
        clear cip
        C = [C; ContourIndPing_i(j,:)];
    end
    
end




clear all
close all
addpath(".\gfunctions")
tic
% Load Label.mat
path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\Label\Label.mat';
load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;

%Remove Removed Contours
ContourDepth(RemovedContour,:) =[];
ContourIndPing(RemovedContour,:) =[];
LabelIndex(RemovedContour) =[];
ContourIndDepth = NaN(size(ContourDepth));
% ContourIndPingEI = NaN(size(ContourDepth));


%load EI features
filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
Classif = matfile(filename,'Writable',false); clear filename;
EsuPingStartEnd = Classif.EsuPingStartEnd;
depth_surface = Classif.depth_surface;
maskFish = Classif.Mask_clean;
Bottom = Classif.Bottom;
DateTime = Classif.DateTime;
boolean_maskFish = ~isnan(maskFish(:,:,1)); %Set boolean matrix echogram ping in x and depth in y axis
Sv_surface = Classif.Sv_surface;
Sv_Fish = Sv_surface.*maskFish;
Sv_Fish = Sv_Fish(:,:,1);


% get CountourDepth index
for i=1:size(ContourDepth,1)
    for j=1:size(ContourDepth,2)
        ind = find(abs(depth_surface-ContourDepth(i,j))==min(abs(depth_surface-ContourDepth(i,j))));
        if ~isempty(ind)
            ContourIndDepth(i,j) = ind(1);
        end
    end
end

% % contour ind ping
% for i=1:size(ContourIndPing,1)
%     for j=1:size(ContourIndPing,2)
% ContourIndPingEI(i,j) = EsuPingStartEnd(1,ContourIndPing(i,j));
%     end
% end

%construct mask

for i=2:max(LabelIndex)
    count = 0 ;
    path_lib = 'C:\Users\jsalveta\Desktop\Matlab code\Labelling\';
    if ~exist([path_lib, LabelName{i}], 'dir')
        mkdir([path_lib, LabelName{i}]);
    end
    ContourIndPing_i = ContourIndPing(LabelIndex==i,:);
    ContourIndDepth_i = ContourIndDepth(LabelIndex==i,:);
    ContourDepth_i = ContourDepth(LabelIndex==i,:);
    
    for j=1:size(ContourIndPing_i,1)
        count = count +1;
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
            f=ContourDepth_i(j,:);
            f=f(~isnan(f));
            f=f(~isnan(cip));
            BW = poly2mask(d,e,size(boolean_maskFish,1),size(boolean_maskFish,2));
            im_bin = boolean_maskFish & BW;
            maskLabel = Sv_Fish.*im_bin;
            gtsecho(maskLabel,depth_surface,[-80 -40],Bottom(1,:));
            hold on
            plot(d,f,'--b','Linewidth',2)
            xlim([min(d)-100 max(d)+100])
            ylim([0 max(max(Bottom(1,min(d)-50:max(d)+50)),max(f))])
            saveas(gcf,[path_lib, LabelName{i},'\FAROFA_2_Label_n',num2str(count),'.png'])
            close all
            clear cip
        else
            clear cip
            C = [C; ContourIndPing_i(j,:)];
        end
        
    end
    toc
    % figure;imagesc(maskLabel)
    % figure;imagesc(boolean_maskFish)
    % figure;imagesc(BW)
    % figure;imagesc(im_bin)
    % figure;imagesc(Sv_Fish(:,:,1))
    % gtsecho(Sv_Fish(:,:,1),depth_surface,[-80 -40],Bottom(1,:));
    %
    
    
end