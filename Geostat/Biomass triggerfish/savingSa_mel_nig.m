clear all
close all
addpath("C:\Users\jsalveta\Desktop\Matlab code\gfunctions")


%for each Farofa
for FAROFA=1:3;
    for is_video=1:2% 2 pou melange 
        if FAROFA == 1
            path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\Label\Label.mat';
            filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_1\EK80\HAC\Cruise_FAROFA_1\Treatment20191128_201456_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';
            if is_video == 0
                label_i=12;
            elseif is_video == 1
                label_i = 17;
            elseif is_video == 2
                label_i = [21, 22];
            end
            
        elseif FAROFA == 2
            path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\Label\Label.mat';
            filename ='C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\CleanResults\Echointegration\EI001\Echointegration.mat';            
            if is_video == 0
                label_i=12;
            elseif is_video == 1
                label_i =4;
            elseif is_video == 2
                label_i = [28,29,50];
            end
            
        elseif FAROFA ==3
            path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_3\Label\Label.mat';
            filename = 'C:\Users\jsalveta\Desktop\Data\FAROFA_3\Echointegration_FISH_cleaned\Echointegration.mat';
            if is_video == 0
                label_i=12;
            elseif is_video == 1
                label_i =17;
            elseif is_video == 2
                label_i = [21,27,28];
            end
        end
        % Load Label.mat
        load(path2Label,'ContourDepth','ContourIndPing','LabelIndex','LabelName','RemovedContour'); clear path2Label;
        
        ContourDepth(RemovedContour,:) = [];
        ContourIndPing(RemovedContour,:) = [];
        LabelIndex(RemovedContour) = [];
        ContourIndDepth = NaN(size(ContourDepth));
        
        
        %load EI features
        Classif = matfile(filename,'Writable',false); clear filename;
        depth_surface = Classif.depth_surface;
        maskFish = Classif.Mask_clean(:,:,1);
        Time = (Classif.Time / (24*3600))+datenum('01-Jan-1970');
        EsuPingStartEnd = Classif.EsuPingStartEnd;
        [Y M D]=datevec(Time);
        Bottom = Classif.Bottom(1,:);
        Latitude = Classif.Latitude;
        Longitude = Classif.Longitude;
        
        SvFish=Classif.Sv_surface(:,:,1).*maskFish;
        Sa_surface_fish_70 = Classif.Sa_surface(:,:,1).*maskFish;
        Sa_surface_fish_200 = Classif.Sa_surface(:,:,2).*maskFish;
        clear maskFish
        maskBottom = zeros(size(Sa_surface_fish_70));
for i=1:length(Bottom)
    [minValue,closestIndex2] = min(abs(depth_surface-(Bottom(i)-1.5)));
    maskBottom(closestIndex2+1:end,i)=ones(size(maskBottom,1)-closestIndex2,1);
end

        
        
        boolean_maskFish = ~isnan(SvFish(:,:,1)) & ~maskBottom; %Set boolean matrix echogram ping in x and depth in y axis
clear maskBottom SvFish
        %         figure; imagesc(boolean_maskFish(1:400,1:10000))
%         figure; imagesc(Sa_surface_fish_70(1:400,1:10000))
%  figure;gtsecho(Sa_surface_fish_70(1:400,1:10000),depth_surface(1:400),[0 40],Bottom(1:10000));title('Original EI 70 kHz');

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
        
        BWf = zeros(size(boolean_maskFish));
        
        for l=1:length(label_i)
        ContourIndPing_i = ContourIndPing(LabelIndex==label_i(l),:);
        ContourIndDepth_i = ContourIndDepth(LabelIndex==label_i(l),:);
        
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
        end
        im_bin = boolean_maskFish & BWf;
        %maskLabel(im_bin)=1;
        clear maskFish maskLabel BW BWf boolean_maskFish SvFish
        %end
        % figure;imagesc(maskLabel)
        % figure;imagesc(im_bin)
        Sa_mel_nig_vid_70=gmask(Sa_surface_fish_70,im_bin,0);
        Sa_mel_nig_vid_200=gmask(Sa_surface_fish_200,im_bin,0);
        Horizontal_sa_surface_fish_70 = squeeze(nansum(Sa_mel_nig_vid_70,1));
        Horizontal_sa_surface_fish_200 = squeeze(nansum(Sa_mel_nig_vid_200,1));
        
         figure;gtsecho(Sa_surface_fish_70(:,57312:59150),depth_surface',[0 40],Bottom(57312:59150));colorbar;title('Sa fish 70 kHz');ylim([0 100])
        % figure; imagesc(Sa_surface_fish_70(:,57312:59150))
        % figure; imagesc(Sa_mel_nig_vid_70)
        
        EsuPingStart=EsuPingStartEnd(1,:);
        %video = ones(size(Latitude));
        video = repmat(is_video,size(Latitude));
        F = repmat(FAROFA,size(Latitude));
        offset_projection = 0.05; %in degree of Lat and long : extend boudaries of map
        
        Longitude_min = -32.65;
        Longitude_max = -32.349;
        Latitude_min = -3.91;
        Latitude_max = -3.79;
        diff_lat=[0 abs(diff(Latitude))>0.001];
        diff_lon=[0 abs(diff(Longitude))>0.001];
        out_boundary=find(Latitude>Latitude_max+offset_projection| Latitude<Latitude_min-offset_projection|Longitude>Longitude_max+offset_projection| Longitude<Longitude_min-offset_projection| diff_lon==1|diff_lat==1);
        %figure;plot(Longitude,Latitude,'b.',Longitude(out_boundary),Latitude(out_boundary),'r.')
        
        Longitude(out_boundary)=[];
        Latitude(out_boundary)=[];
        Horizontal_sa_surface_fish_70(out_boundary)=[];
        Horizontal_sa_surface_fish_200(out_boundary)=[];
        Time(out_boundary)=[];
        Bottom(out_boundary)=[];
        EsuPingStart(out_boundary)=[];
        video(out_boundary)=[];
        F(out_boundary)=[];
        
        T = table(Longitude',Latitude',Horizontal_sa_surface_fish_70',Horizontal_sa_surface_fish_200',Time',Bottom',EsuPingStart',F',video','VariableNames',{'Longitude','Latitude','Horizontal_sa_surface_fish_70','Horizontal_sa_surface_fish_200','Time','Bottom','EsuPingStart','FAROFA','Video'});
        writetable(T,['F',num2str(FAROFA),'_mel_nig_vid_',num2str(is_video),'.csv'])
    end
end
