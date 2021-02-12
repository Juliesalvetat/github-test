%%%%%%% Calculate Sa horizontal vertical sum for the different EI (fish, no
%%%%%%% fish, original in each folder

clear all
addpath("C:\Users\jsalveta\Desktop\Matlab code\gfunctions")

folder = {'FISH' 'NO FISH' 'ORIGINAL'};

for f=1:length(folder)
    
    Class_Filename = ['.\',folder{f},'\Echointegration.mat'];
    corr_Path = ['.\',folder{f},'\Sa_Horizontal.mat'];
    
    Classif = matfile(Class_Filename,'Writable',false);
    corr = matfile(corr_Path,'Writable',true);
    
    Esu_step = 1000;
    totalNesu=size(Classif.depth_bottom,2);
    seq_ESU=1:Esu_step:totalNesu;
    if  seq_ESU(end)== totalNesu
        seq_ESU=seq_ESU;
    else
        seq_ESU=[seq_ESU,totalNesu];
    end
    for i= 1:length(seq_ESU)-1
        disp([num2str(i),'/',num2str(length(seq_ESU)),' ',folder{f}])
        
        if i==1
            %Mask_clean = Classif.Mask_clean(:,seq_ESU(i):seq_ESU(i+1),:);
            corr.Sa_surface = Classif.Sa_surface(:,seq_ESU(i):seq_ESU(i+1),:);%.*Mask_clean;
            corr.Sv_surface = Classif.Sv_surface(:,seq_ESU(i):seq_ESU(i+1),:);%.*Mask_clean;
            
            corr.Horizontal_sa_surface = squeeze(nansum(corr.Sa_surface,1));
            corr.Horizontal_sv_surface = 10.*log10(squeeze(nanmean(10.^(corr.Sv_surface),1)));
            
        else
            Mask_clean = Classif.Mask_clean(:,seq_ESU(i)+1:seq_ESU(i+1),:);
            corr.Sa_surface(:,seq_ESU(i)+1:seq_ESU(i+1),:) = Classif.Sa_surface(:,seq_ESU(i)+1:seq_ESU(i+1),:);%.*Mask_clean;
            corr.Sv_surface(:,seq_ESU(i)+1:seq_ESU(i+1),:) = Classif.Sv_surface(:,seq_ESU(i)+1:seq_ESU(i+1),:);%.*Mask_clean;
            
            corr.Horizontal_sa_surface(seq_ESU(i)+1:seq_ESU(i+1),:) = squeeze(nansum(corr.Sa_surface(:,seq_ESU(i)+1:seq_ESU(i+1),:),1));
            corr.Horizontal_sv_surface(seq_ESU(i)+1:seq_ESU(i+1),:) = 10.*log10(squeeze(nanmean(10.^(corr.Sv_surface(:,seq_ESU(i)+1:seq_ESU(i+1),:)),1)));
            %             end
        end
        
        
    end
end
    % figure;imagesc(Classif.Sa_surface(:,seq_ESU(i):seq_ESU(i+1),1))
    % SaSum = Classif.SaSum;
    % Horizontal_sa_surface = corr.Horizontal_sa_surface;
    % addpath("C:\Users\jsalveta\Desktop\Matlab code\gfunctions")
    % gtsecho(Classif.Sa_surface(:,seq_ESU(i):seq_ESU(i+1),1),Classif.depth_surface,[0 40],Classif.Bottom(1,seq_ESU(i):seq_ESU(i+1)));title('Original EI 70 kHz');
    % ylim([0 max(Classif.Bottom(1,seq_ESU(i):seq_ESU(i+1)))+1])
    %
    %
    %
    % Longitude = Classif.Longitude;
    % Latitude = Classif.Latitude;
    % DateTime = Classif.DateTime;
    %
    % figure;plot(Longitude,Latitude,'.')
    %
    % offset_projection = 0.05; %in degree of Lat and long : extend boudaries of map
    % Longitude_min = -32.482;
    % Longitude_max = -32.349;
    % Latitude_min = -3.91;
    % Latitude_max = -3.79;
    %
    % addpath("C:\Users\jsalveta\Desktop\Matlab code\m_map\")
    % figure
    % m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    % m_pcolor(Xgrid,Ygrid,squeeze(log10(distrib_sa_surface_fish(:,:,1)+1))');
    % shading flat;
    % m_gshhs_f('patch',[.5 .5 .5]);
    % m_grid('linestyle','none','box','fancy','fontsize',14,'fontweight','bold');
    %
    %
    % diff_lat=[0 abs(diff(Latitude))>0.01];
    % diff_lon=[0 abs(diff(Longitude))>0.01];
    % out_boundary=find(Latitude>Latitude_max+offset_projection| Latitude<Latitude_min-offset_projection|Longitude>Longitude_max+offset_projection| Longitude<Longitude_min-offset_projection| diff_lon==1|diff_lat==1);
    %
    % Longitude(out_boundary)=[];
    % Latitude(out_boundary)=[];
    % figure;plot(Longitude,Latitude,'b.',Longitude(out_boundary),Latitude(out_boundary),'r.')
