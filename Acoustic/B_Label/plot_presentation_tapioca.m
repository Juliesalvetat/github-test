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

load('C:\Users\jsalveta\Desktop\Matlab code\bathy_coast_FdN\poubelle\Bati_Ilhas.mat')
[C,h]=m_contour(lon,lat,Elevation,[-80,-3000],'linewi',1,'color','k')%,'showtext','on','labelspacing',3000);
  m_line(X1,X2,'color','r','linewi',1.5) % MPA limits
         m_hatch(X1,X2,'single',30,5,'color','r'); % ...with hatching added.
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


figure
    m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    m_usercoast('FN','patch',[.5 .5 .5]); hold on % can't remember how I did that savec FN coastline
    % must something like
    % %m_gshhs_f('patch',[.5 .5 .5]);
    % %m_gshhs_h('save','FN');
    % m_grid('box','fancy','fontsize',14,'fontweight','bold');
    % hold on
     m_line(X1,X2,'color','r','linewi',1.5) % MPA limits
         m_hatch(X1,X2,'single',30,5,'color','r'); 
repres_line_y = [y_lineW(x_line<=X(2)) y_lineE(x_line>=X(2))];
        m_line(x_line,repres_line_y,'color','b','linewi',1.5)
           h1=m_line(Poly_lon(end),Poly_lat(end),'marker','o','color',[0 .5 0],...
          'linest','none','markerfacecolor','r','clip','point');
    h1=m_line(Poly_lon(end),Poly_lat(end),'marker','o','color',[.49 1 .63],...
          'linest','none','LineWidth',2,...
                                                     'MarkerEdgeColor','k',...
                                                     'MarkerFaceColor',[.49 1 .63],...
                                                     'MarkerSize',10)
load('FN')
long = ncst(:,1);
lat = ncst(:,2);
coastLong = long(~isnan(long));
    coastLat = lat(~isnan(lat));
    
    lat = Poly_lon(end);        %# Inland latitude points (in degrees)
    long = Poly_lat(end);  %# Inland longitude points (in degrees)
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
    
    
    plot(coastLong,coastLat,'.');hold on
    for i=1:length(long)
        
  plot([m_ll2xy(long(i))  m_ll2xy(coastLong(coastIndex(i)))],...    %'# Plot the inland points and
         [m_ll2xy(lat(i)) m_ll2xy(coastLat(coastIndex(i)))],...      %'#   nearest coastal points
         'ro-');
    end
        h1=m_line(coastLong(coastIndex(i)),coastLat(coastIndex(i)),'marker','o','color',[0 .5 0],...
          'linest','none','markerfacecolor','r','clip','point');