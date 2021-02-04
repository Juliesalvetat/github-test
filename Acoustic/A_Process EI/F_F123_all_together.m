clear all; close all
addpath("C:\Users\jsalveta\Desktop\Matlab code\m_map\")
offset_projection = 0.05; %in degree of Lat and long : extend boudaries of map
Longitude_min = -32.482;
Longitude_max = -32.349;
Latitude_min = -3.91;
Latitude_max = -3.79;

% Longitude_min = -33.5;
% Longitude_max = -31.5;
% Latitude_min = -4.6;
% Latitude_max = -3;


folder = {'FISH' 'NO FISH' 'ORIGINAL'};

for fo=1:length(folder)
    LAT = [];
LON = [];
DEPTH = [];
H_fish = [];
H_tot = [];
TIME = [];
F = [];
ESU = [];
    for FAROFA=1:3
disp(['FAROFA ',num2str(FAROFA),' : ',folder{fo}])
    
        Class_Filename = ['.\FAROFA',num2str(FAROFA),'\',folder{fo},'\Echointegration.mat'];
        corr_Path = ['.\FAROFA',num2str(FAROFA),'\',folder{fo},'\Sa_Horizontal.mat'];
        
        Classif = matfile(Class_Filename,'Writable',false);
        corr = matfile(corr_Path,'Writable',true);
        Horizontal_sa_surface = corr.Horizontal_sa_surface;
        
        depth_bottom = Classif.depth_bottom(:,:,1);
        Longitude=Classif.Longitude;
        Latitude=Classif.Latitude;
        Time = (Classif.Time / (24*3600))+datenum('01-Jan-1970');
        f = repmat(FAROFA,size(Time));
        c=Classif.EsuPingStartEnd(1,:);        
        LAT = [LAT Latitude];
        LON = [LON Longitude];
        DEPTH = [DEPTH depth_bottom];
        H_tot = [H_tot; Horizontal_sa_surface];
        TIME = [TIME Time];
        F = [F f];
        ESU = [ESU c];
        
        
        %
        % Horizontal_sv_surface200 =Horizontal_sv_surface(~isnan(Horizontal_sa_surface(:,2)),2);
        % Horizontal_sa_surface200 =Horizontal_sa_surface(~isnan(Horizontal_sa_surface(:,2)),2);
        % Horizontal_log_sa_surface200 =Horizontal_log_sa_surface(~isnan(Horizontal_sa_surface(:,2)),2);
        % Longitude=Longitude(~isnan(Horizontal_sa_surface(:,2)));
        % Latitude=Latitude(~isnan(Horizontal_sa_surface(:,2)));
        % depth_bottom = depth_bottom(~isnan(Horizontal_sa_surface(:,2)));
        %
        % figure;scatter(Longitude,Latitude,[],Horizontal_sv_surface200,'filled')
        % figure;scatter(Longitude,Latitude,[],Horizontal_sa_surface200,'filled')
        % figure;scatter(Longitude,Latitude,[],Horizontal_log_sa_surface200,'filled')
        
        %    T = table(Longitude',Latitude',depth_bottom',Horizontal_sv_surface200,Horizontal_sa_surface200,Horizontal_log_sa_surface200,'VariableNames',{'Longitude','Latitude','depth_bottom','Horizontal_Sv200','Horizontal_Sa200','Horizontal_logSa200'});
        % writetable(T,'FAROFA1_Horizontal_Sa.csv')
        
        
        

        
        % figure;plot(Longitude,Latitude,'.')
    end
        Latitude = LAT;
        Longitude =  LON;
        Horizontal_sa_surface = H_tot;
        
        Depth = DEPTH;
        Time = TIME;
        
        
        diff_lat=[0 abs(diff(Latitude))>0.01];
        diff_lon=[0 abs(diff(Longitude))>0.01];
        out_boundary=find(Latitude>Latitude_max+offset_projection| Latitude<Latitude_min-offset_projection|Longitude>Longitude_max+offset_projection| Longitude<Longitude_min-offset_projection| diff_lon==1|diff_lat==1);
        
        if (~isempty(out_boundary))
            Latitude(out_boundary)=[];
            Longitude(out_boundary)=[];
            Depth(out_boundary)=[];
            Horizontal_sa_surface(out_boundary,:)=[];
            
            Time(out_boundary)=[];
            F(out_boundary)=[];
            ESU(out_boundary)=[];
        end
    save(['./FAROFA123/Sa_Horizontal_F123_',folder{fo},'.mat'],'Longitude','Latitude','Depth','Horizontal_sa_surface','Time','F','ESU')
     T = table(Longitude',Latitude',Horizontal_sa_surface(:,1),Horizontal_sa_surface(:,2),Depth',Time',ESU',F','VariableNames',{'Longitude','Latitude','Horizontal_Sa70','Horizontal_Sa200','Bottom_Depth','Time','ESU','FAROFA'});
 writetable(T,['./FAROFA123/Sa_Horizontal_F123_1ping_',folder{fo},'.csv'])
end
datestr(Time(Horizontal_sa_surface(:,1)>100000))
figure;plot(Longitude,Latitude,'.')
hold on
plot(Longitude(Horizontal_sa_surface(:,1)>10000),Latitude(Horizontal_sa_surface(:,1)>10000),'.r')
%
% %ABRACOS 2
% Class_Filename = 'C:\Users\jsalveta\Desktop\Data\ABRACOS 2\LEG 2\Echointegration.mat';
% corr_Path = '.\Sa_Horizontal_A2.mat';
%
% Horizontal_sa_surfaceA = corr.Horizontal_sa_surface;
%
% % Horizontal_sv_surface_fish = corr.Horizontal_sv_surface_fish;
%
% depth_bottomA = Classif.depth_bottom(:,:,1);
% LongitudeA=Classif.Longitude;
% LatitudeA=Classif.Latitude;
% TimeA = (Classif.Time / (24*3600))+datenum('01-Jan-1970');
% f = repmat(FAROFA,size(Time));
% c=1:length(TimeA);
%
%
% Night1Sunrise2Day3Sunset4=Classif.Night1Sunrise2Day3Sunset4;
% LatitudeA = LatitudeA(Night1Sunrise2Day3Sunset4==3);
% LongitudeA = LongitudeA(Night1Sunrise2Day3Sunset4==3);
%
% % somme sur la dimension vertical
%
% offset_projection = 0.05; %in degree of Lat and long : extend boudaries of map
% Longitude_min = -33.3;
% Longitude_max = -31.55;
% Latitude_min = -4.6;
% Latitude_max = -3;
%
% diff_lat=[0 abs(diff(Latitude))>0.1];
% diff_lon=[0 abs(diff(Longitude))>0.1];
% out_boundary=find(Latitude>Latitude_max+offset_projection| Latitude<Latitude_min-offset_projection|Longitude>Longitude_max+offset_projection| Longitude<Longitude_min-offset_projection| diff_lon==1|diff_lat==1);
%
% writetable(T,'FAROFA123_logSa_1ping.csv')
%
% T=readtable('FAROFA123_logSa_1ping.csv');
%
vect_grid_size=  [50 100 250 500 1000];
%diff = Horizontal_sa_surface - Horizontal_sa_surface_fish;
% diff(isnan(diff))=-Inf;
% [M,I]=sort(diff,'descend');
% figure;plot(Longitude,Latitude,'.');hold on
% plot(Longitude(I(1:100,1)),Latitude(I(1:100,1)),'.r')
% figure;boxplot(diff)

g=3;
%for g=2:length(vect_grid_size)-1
    grid_size = vect_grid_size(g)/(1852*60);
    dlon = grid_size;
    dlat = grid_size;
    Xgrid = Longitude_min-offset_projection:dlon:Longitude_max+offset_projection;
    Ygrid = Latitude_min-offset_projection:dlat:Latitude_max+offset_projection;
    [xx,yy] = meshgrid(Xgrid,Ygrid);
    
    
    distrib_sa_surface_fish = NaN(length(Xgrid),length(Ygrid),2);
    distrib_sa_surface = NaN(length(Xgrid),length(Ygrid),2);
    
    distrib_sa_diff = NaN(length(Xgrid),length(Ygrid),2);
    
    Bottom = NaN(length(Xgrid),length(Ygrid),2);
    Lat = NaN(length(Xgrid),length(Ygrid),2);
    Lon = NaN(length(Xgrid),length(Ygrid),2);
    
    for indice_frequence=1:2
        count=0;
        for i=1:length(Xgrid)
            for j=1:length(Ygrid)
                count=count+1;
                disp([num2str(count),'/',num2str(length(Xgrid)*length(Ygrid)),' freq',num2str(indice_frequence)])
                %indices d'ESU correspondant au point de grille
                ind = find(Longitude>Xgrid(i) & Longitude<=Xgrid(i)+dlon & Latitude>Ygrid(j) & Latitude<=Ygrid(j)+dlat);
                if (~isempty(ind))
                    distrib_sa_surface_fish(i,j,indice_frequence) = nanmean(Horizontal_sa_surface_fish(ind,indice_frequence),1);
                    distrib_sa_surface(i,j,indice_frequence) = nanmean(Horizontal_sa_surface(ind,indice_frequence),1);
                    distrib_sa_diff(i,j,indice_frequence) = nanmean(diff(ind,indice_frequence),1);
                    
                    Bottom(i,j,indice_frequence) = nanmean(Depth(ind));
                    Lat(i,j,indice_frequence) = nanmean(Latitude(ind));
                    Lon(i,j,indice_frequence) = nanmean(Longitude(ind));
                end
            end
        end
    end
    
    
    %     diff =  distrib_sa_surface - distrib_sa_surface_fish;
    
    Sa_2D_70=distrib_sa_surface(:,:,1);
    Sa_2D_200=distrib_sa_surface(:,:,2);
    Sa_2D_70_FISH=distrib_sa_surface_fish(:,:,1);
    Sa_2D_200_FISH=distrib_sa_surface_fish(:,:,2);
    Sa_2D_70_NO_FISH=distrib_sa_diff(:,:,1);
    Sa_2D_200_NO_FISH=distrib_sa_diff(:,:,2);
    Sa_2D_70=Sa_2D_70(:);
    Sa_2D_200=Sa_2D_200(:);
    Sa_2D_70_FISH=Sa_2D_70_FISH(:);
    Sa_2D_200_FISH=Sa_2D_200_FISH(:);
    Sa_2D_70_NO_FISH=Sa_2D_70_NO_FISH(:);
    Sa_2D_200_NO_FISH=Sa_2D_200_NO_FISH(:);
    
    LON_Sa_2D = Lon(:,:,1);
    LAT_Sa_2D = Lat(:,:,1);
    LON_Sa_2D = LON_Sa_2D(:);
    LAT_Sa_2D = LAT_Sa_2D(:);
    
    Bottom = Bottom(:,:,1);
    Bottom = Bottom(:);
    
    % figure;plot(LON_Sa_2D,LAT_Sa_2D,'.')
    % ESU = 1:length(Bottom);
    T = table(LON_Sa_2D,LAT_Sa_2D,Sa_2D_70,Sa_2D_200,Sa_2D_70_FISH,Sa_2D_200_FISH,Sa_2D_70_NO_FISH,Sa_2D_200_NO_FISH,Bottom,'VariableNames',{'Longitude','Latitude','Horizontal_Sa70_total','Horizontal_Sa200_total','Horizontal_Sa70_fish','Horizontal_Sa200_fish','Horizontal_Sa70_no_fish','Horizontal_Sa200_no_fish','Bottom_Depth'});
    writetable(T,'FAROFA123_logSa_250m_fuller.csv')
    figure
    %set(gcf, 'Position', get(0, 'Screensize'));
    % subplot(1,2,1)
    m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    m_pcolor(Xgrid,Ygrid,squeeze(log10(distrib_sa_surface_fish(:,:,1)+1))');
    %m_pcolor(Xgrid,Ygrid,squeeze(log10(distrib_sa_diff(:,:,1)+1))');
    
    shading flat;
    m_gshhs_f('patch',[.5 .5 .5]);
    m_grid('linestyle','none','box','fancy','fontsize',14,'fontweight','bold');
    h=colorbar;
    caxis([0 2])
    m_line(X1,X2,'linewi',2,'color','r');
    
    log_distrib_log_sa_surface_cell_fish{g,1} = log_distrib_log_sa_surface_fish;
    log_distrib_log_sa_surface_cell_fish{g,2} = vect_grid_size(g);
    
    
    %  corr.log_distrib_log_sa_surface = log_distrib_log_sa_surface;
    %  corr.distrib_leg_sa_surface = distrib_leg_sa_surface;
    %  corr.distrib_leg_sv_surface = distrib_leg_sv_surface;
    %
    %   corr.log_distrib_log_sa_surface_fish = log_distrib_log_sa_surface_fish;
    %  corr.distrib_leg_sa_surface_fish = distrib_leg_sa_surface_fish;
    %  corr.distrib_leg_sv_surface_fish = distrib_leg_sv_surface_fish;
end
g = 2;
grid_size = vect_grid_size(g)/(1852*60);
dlon = grid_size;
dlat = grid_size;
Xgrid = Longitude_min-offset_projection:dlon:Longitude_max+offset_projection;
Ygrid = Latitude_min-offset_projection:dlat:Latitude_max+offset_projection;
[xx,yy] = meshgrid(Xgrid,Ygrid);
temp = log_distrib_log_sa_surface_cell_fish{g,1};
t = temp(:,:,2);
t(isnan(t))=[];
Vq = interp2(Xgrid,Ygrid,temp(:,:,2)',xx',yy');
Vq = interp2(temp(:,:,2));
Vq = interp2(temp(:,:,2)');

a = isnan(temp(:,:,2));
X=xx(~a);
Y=yy(~a);
t = temp(:,:,2);
t = t(~a);
Vq = interp2(X,Y,t);,xx(:)',yy(:)');
Vq = interp2(temp(:,:,2));

figure
set(gcf, 'Position', get(0, 'Screensize'));
% subplot(1,2,1)
m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
m_pcolor(Xgrid,Ygrid,Vq');
shading flat;
m_gshhs_f('patch',[.5 .5 .5]);
m_grid('linestyle','none','box','fancy','fontsize',14,'fontweight','bold');
h=colorbar;

figure;imagesc(Vq)
temp = log_distrib_log_sa_surface_cell_fish{g,1};
temp=temp(:,:,2);
t = ~isnan(temp);
[m,n] = size(temp);
F=scatteredInterpolant(xx(t),yy(t),temp(t))
[X,Y]=ndgrid(0:0.1:m-1,0:0.1:n-1);
figure; mesh(X,Y,F(X,Y))

corr.log_distrib_log_sa_surface_cell = log_distrib_log_sa_surface_cell;
corr.distrib_leg_sa_surface_cell = distrib_leg_sa_surface_cell;
corr.distrib_leg_sv_surface_cell = distrib_leg_sv_surface_cell;

corr.log_distrib_log_sa_surface_cell_fish = log_distrib_log_sa_surface_cell_fish;
corr.distrib_leg_sa_surface_cell_fish = distrib_leg_sa_surface_cell_fish;
corr.distrib_leg_sv_surface_cell_fish = distrib_leg_sv_surface_cell_fish;
end