clear all;close all;clc

%%%%% Add some stuff to Polygone_multivar table
% MPA
% wind exposition
% shelfbreak distance
% coast distance

% add m_map path
addpath('C:\Users\jsalveta\Desktop\Matlab code\m_map\')
% set boundarie limits
%Longitude_min = -32.65; (with Drina submount)
Longitude_min = -32.5; % (without Drina submount)
Longitude_max = -32.349;
Latitude_min = -3.91;
Latitude_max = -3.79;
offset_projection = 0.005;
load('cursor_info_new.mat')

% load that PARNAMAR_limits.mat containing PARNAMAR limits
load('PARNAMAR_limits.mat')
figure
m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
m_usercoast('FN','patch',[.5 .5 .5]); hold on
m_grid('linestyle','none','fontsize',14,'fontweight','bold');hold on

cursor_info = cursor_info_new;
x1=cursor_info(1).Position(1);x2=cursor_info(2).Position(1);x3=cursor_info(3).Position(1);
y1=cursor_info(1).Position(2);y2=cursor_info(2).Position(2);y3=cursor_info(3).Position(2);
line([x1 x2],[y1 y2],'color','r','linewi',1.5)
line([x2 x3],[y2 y3],'color','r','linewi',1.5)
x = [cursor_info(1).Position(1) cursor_info(2).Position(1) cursor_info(3).Position(1)];
y = [cursor_info(1).Position(2) cursor_info(2).Position(2) cursor_info(3).Position(2)];
[X,Y] = m_xy2ll(x,y);
coefficients_west = polyfit(X(2:3),Y(2:3), 1);
aW = coefficients_west (1);
bW = coefficients_west (2);
coefficients_east = polyfit(X(1:2),Y(1:2), 1);
aE = coefficients_east (1);
bE = coefficients_east (2);
x_line = [Longitude_min:0.001:Longitude_max];
y_lineW = aW*x_line+bW;
y_lineE = aE*x_line+bE;
m_line(
load('FN')
long = ncst(:,1);
lat = ncst(:,2);


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


for FAROFA=2:3
    % read Polygone table from Polygone_multivar.m
    Poly_table = readtable(['F',num2str(FAROFA),'_masklabel.csv']);
    Longitude = Poly_table.Poly_Longitude;
    Latitude = Poly_table.Poly_Latitude;
    Bottom = Poly_table.Poly_Bottom;
    
    %%%%%%%%%%%%%%%%%%%% to be inside the MPA or not to be (inside the MPA) ?
    %%%%%% how do you do ? (glad you asked)
    %%% there is this file PARNAMAR_limits.mat that containes the MPA limits
    %%% it was made from a klm of the PARNAMAR limits (thx Sophie :D)
    %%% not sure how I did the transformation from klm to PARNAMAR_limits.mat
    %%% file but it work so who cares? do you ? -_-" ok I'll find it somewhere
    %%% somehow
    
    %%% so basicaly it's just testing if lat and lon are inside or outside the
    %%% MPA BUT at equals bottom depths meaning <50m (plot twist)
    
    %%% wana see ?
    
    %%% try this plot
    
    figure
    m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    m_usercoast('FN','patch',[.5 .5 .5]); hold on % can't remember how I did that savec FN coastline
    % must something like
    % %m_gshhs_f('patch',[.5 .5 .5]);
    % %m_gshhs_h('save','FN');
    % m_grid('box','fancy','fontsize',14,'fontweight','bold');
    % hold on
     m_line(X1,X2,'color','k','linewi',1.5) % MPA limits
         m_hatch(X1,X2,'single',30,5,'color','k'); % ...with hatching added.
    % m_line(Longitude,Latitude,'marker','.','color',[.5 .5 .5],'linewi',1,...
    %           'linest','none','markerfacecolor','[.5 .5 .5]');
    
    %%% so now you wanna know what data is above 50m and inside or outside the
    %%% PARNAMAR
    
    in = inpolygon(Longitude,Latitude,X1,X2); % data in or on polygone
    
    % quick visual check
    % m_line(Longitude(in),Latitude(in),'marker','.','color','r','linewi',1,...
    %            'linest','none','markerfacecolor','r')
    % m_line(Longitude(~in & Bottom<=50),Latitude(~in & Bottom<=50),'marker','.','color','b','linewi',1,...
    %            'linest','none','markerfacecolor','b')
    % m_line(Longitude(~in & Bottom>50),Latitude(~in & Bottom>50),'marker','.','color','g','linewi',1,...
    %            'linest','none','markerfacecolor','g')
    %
    %
    % now you want to add that information to your table, but how?
    % glad you asked
    % lets have a code 'in' : for inside MPA; 'out' : for outside MPA; 'off50m'
    % : for outside 50m depth
    
    MPA = strings(length(Longitude),1);
    MPA(in)='in';
    MPA(~in & Bottom<=50) = 'out';
    MPA(~in & Bottom>50) = 'off50m';
    
    % done :)
    
    %%% now you want to know if you are on the exposed side of the island or
    %%% not be be leeward or to be windward?
    
    % first you need to define what is windward what is leeward ?
    % y drawing a line on figure Tools> DataCursor click and push shift to
    % create new points then right click select export Cursor data to workspace
    % save cursor_info_new.mat cursor_info_new
    % cursor_info = cursor_info_new;
    % x1=cursor_info(1).Position(1);x2=cursor_info(2).Position(1);x3=cursor_info(3).Position(1);
    % y1=cursor_info(1).Position(2);y2=cursor_info(2).Position(2);y3=cursor_info(3).Position(2);
    % % line([x1 x2],[y1 y2],'color','r','linewi',1.5)
    % % line([x2 x3],[y2 y3],'color','r','linewi',1.5)
    % x = [cursor_info(1).Position(1) cursor_info(2).Position(1) cursor_info(3).Position(1)];
    % y = [cursor_info(1).Position(2) cursor_info(2).Position(2) cursor_info(3).Position(2)];
    % [X,Y] = m_xy2ll(x,y);
    % m_line(X,Y,'color','k','linewi',1.5)
    % coefficients_west = polyfit(X(2:3),Y(2:3), 1);
    % aW = coefficients_west (1);
    % bW = coefficients_west (2);
    % coefficients_east = polyfit(X(1:2),Y(1:2), 1);
    % aE = coefficients_east (1);
    % bE = coefficients_east (2);
    % x_line = [Longitude_min:0.001:Longitude_max];
    % y_lineW = aW*x_line+bW;
    % y_lineE = aE*x_line+bE;
    % m_line(x_line,y_lineW,'color','g','linewi',1.5)
    % m_line(x_line,y_lineE,'color','r','linewi',1.5)
    %
    repres_line_y = [y_lineW(x_line<=X(2)) y_lineE(x_line>=X(2))];
        m_line(x_line,repres_line_y,'color','b','linewi',1.5)
    h1=m_line(Poly_lon(end),Poly_lat(end),'marker','s','color',[0 .5 0],...
          'linest','none','markerfacecolor','r','clip','point');
    %%% we use 2 line portions and their intersection at X(2)
    ind_SousLeVent = find((Longitude<X(2) & Latitude-(aW*Longitude+bW) >= 0) ...
        | (Longitude>X(2) & Latitude-(aE*Longitude+bE) >= 0));
    
    ind_AuVent = find((Longitude<X(2) & Latitude-(aW*Longitude+bW)<0) ...
        |  (Longitude>X(2) & Latitude-(aE*Longitude+bE)<0));
    
    % %%%visual check
    % m_line(Longitude(ind_SousLeVent),Latitude(ind_SousLeVent),'marker','.','color','g','linewi',1,...
    %            'linest','none','markerfacecolor','g')
    % m_line(Longitude(ind_AuVent),Latitude(ind_AuVent),'marker','.','color','b','linewi',1,...
    %            'linest','none','markerfacecolor','b')
    
    %%% add to table
    
    Wind_exposure = strings(length(Longitude),1);
    Wind_exposure(ind_AuVent)='windward';
    Wind_exposure(ind_SousLeVent)='leeward';
    
    %%%%%%%%%%% distance a la cote
    %%% on va utiliser dist_from_coast.m
    %%% small twist coast is actually FN
    
    
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
    % figure
    % plot(coastLong,coastLat,'.');hold on
    % for i=1:length(long)
    % plot([long(i)  coastLong(coastIndex(i))],...    %'# Plot the inland points and
    %      [lat(i) coastLat(coastIndex(i))],...      %'#   nearest coastal points
    %      'ro-');
    % end
    
    
    
    
    % str = strcat(num2str(distanceToCoast(i),...  %'# Make text for the distances
    %                      '%0.1f'),{' km'});
    % text(long(i),lat(i),str(i),'Color','r','VerticalAlignment','bottom');  %# Plot the text
    
    
    
    %%%%%%%%%%%% distance au talu 80m
    % load('C:\Users\jsalveta\Desktop\Matlab code\bathy_coast_FdN\poubelle\Bati_Ilhas.mat')
    % Longitude_min = -32.5; % (without Drina submount)
    % Longitude_max = -32.349;
    % Latitude_min = -3.91;
    % Latitude_max = -3.79;
    % offset_projection = 0.05;
    % ind_lon = find(lon > Longitude_min - offset_projection & lon < Longitude_max + offset_projection);
    % ind_lat = find(lat > Latitude_min - offset_projection & lat < Latitude_max + offset_projection);
    % Elevation =Elevation(ind_lat,ind_lon);
    
    
    % %figure
     [C,h] = contour(lon(ind_lon),lat(ind_lat), Elevation,[-80 -80]); %hold on
    % Lvls = h.LevelList;
    %     idxc = find(C(1,:) == Lvls);
    %     Llen = C(2,idxc);
    %     conturc{1} = C(:,idxc(1)+1 : idxc(1)+1+Llen(1)-1);
    %     conturc{2} = C(:,idxc(2)+1 : idxc(2)+1+Llen(2)-1);
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
%     figure
%     plot(coastLong,coastLat,'.');hold on
%     for i=1:length(long)
%         plot([long(i)  coastLong(coastIndex(i))],...    %'# Plot the inland points and
%             [lat(i) coastLat(coastIndex(i))],...      %'#   nearest coastal points
%             'ro-');
%     end
    
    Poly_table.MPA = MPA;
    Poly_table.Wind_exposure = Wind_exposure;
    Poly_table.distanceToCoast = distanceToCoast;
    Poly_table.distanceToShelfbreak = distanceToShelfbreak;
    writetable(Poly_table,['F',num2str(FAROFA),'_masklabel_suite.csv'])
    
end
