clear all
close all
addpath("C:\Users\jsalveta\Desktop\Matlab code\gfunctions")
tic
% Load Label.mat
path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\Label\Label.mat';
load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;

%Remove Removed Contours
ContourDepth(RemovedContour,:) =[];
ContourIndPing(RemovedContour,:) =[];
LabelIndex(RemovedContour) =[];
ContourIndDepth = NaN(size(ContourDepth));
% ContourIndPingEI = NaN(size(ContourDepth));


%load EI features
filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
Classif = matfile(filename,'Writable',false); clear filename;
EsuPingStartEnd = Classif.EsuPingStartEnd;
depth_surface = Classif.depth_surface;
maskFish = Classif.Mask_clean;
Bottom = Classif.Bottom;
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


C=[];
maskLabel = zeros(size(boolean_maskFish));
%construct mask
    path_lib = 'C:\Users\jsalveta\Desktop\Matlab code\Labelling Code\output\FAROFA1\';

for i=1:19%max(LabelIndex)
    count = 0 ;
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
            saveas(gcf,[path_lib, LabelName{i},'\FAROFA_2_Label_',num2str(i),'n',num2str(count),'.png'])
            close all
            clear cip
        else
            clear cip
            C = [C; ContourIndPing_i(j,:)];
        end
        
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


