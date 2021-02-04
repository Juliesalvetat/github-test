clear all
close all
clc;

% Load Label.mat 
path2Label = 'C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200120_010220_FISH\Label\Label.mat';
load(path2Label,'ContourDepth','ContourIndPing','LabelIndex'); clear path2Label;
save(path2Label,'LabelName','-append') 
ContourIndDepth = NaN(size(ContourDepth));
%load EI features 
filename ='D:\ordi Julie\test FAROFA 3\Cruise_CruiseName\Treatment20190927_120900\CleanResults\Echointegration\EI001\Echointegration_FISH\Echointegration.mat';
load(filename,'depth_surface')
% get CountourDepth index 
for i=1:size(ContourDepth,1)
    for j=1:size(ContourDepth,2)
        ind = find(abs(depth_surface-ContourDepth(i,j))==min(abs(depth_surface-ContourDepth(i,j))));
        if ~isempty(ind)
        ContourIndDepth(i,j) = ind(1);
        end
    end
end


% Load TS.mat 
path2TS = 'D:\ordi Julie\test FAROFA 3\Cruise_CruiseName\Treatment20190927_121657_FISH\CleanResults\TS\TS001\TS.mat';
load(path2TS,'TS70','TS200'); clear path2TS;

% 70 kHz
trackLabel70=TS70(22,:);
TSDepth70 = TS70(2,:); %2-Depth from surface (m)
TSDepth70(trackLabel70==0)=[]; %eliminate single TS (Corrected ping track label (0 for single TS))
for i=1:length(TSDepth70)
        ind = find(abs(depth_surface-TSDepth70(i))==min(abs(depth_surface-TSDepth70(i))));
        if ~isempty(ind)
        TSIndDepth70(i) = ind(1);
        end
    end

TSIndPing70 = TS70(1,:);% 1-Ping index (including removed pings)
TSIndPing70 (trackLabel70==0)=[];
TSLabelIndex70 = TS70(26,:);% 26-User species label (0 for TS with no species label)

TScomp70=TS70(4,:); %compensated TS
TScomp70(trackLabel70==0)=[]; %eliminate single TS (Corrected ping track label (0 for single TS))
trackLabel70(trackLabel70==0)=[];
uniquetrackLabel70 = unique(trackLabel70);


% 200 kHz
trackLabel200=TS200(22,:);
TSDepth200 = TS200(2,:); %2-Depth from surface (m)
TSDepth200(trackLabel200==0)=[]; %eliminate single TS (Corrected ping track label (0 for single TS))
for i=1:length(TSDepth200)
        ind = find(abs(depth_surface-TSDepth200(i))==min(abs(depth_surface-TSDepth200(i))));
        if ~isempty(ind)
        TSIndDepth200(i) = ind(1);
        end
    end

TSIndPing200 = TS200(1,:);% 1-Ping index (including removed pings)
TSIndPing200 (trackLabel200==0)=[];
TSLabelIndex200 = TS200(26,:);% 26-User species label (0 for TS with no species label)

TScomp200=TS200(4,:); %compensated TS
TScomp200(trackLabel200==0)=[]; %eliminate single TS (Corrected ping track label (0 for single TS))
trackLabel200(trackLabel200==0)=[];
uniquetrackLabel200 = unique(trackLabel200);


figure;
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colros.

for i=1:max(LabelIndex)

ContourIndPing_i = ContourIndPing(LabelIndex==i,:);
ContourIndDepth_i = ContourIndDepth(LabelIndex==i,:);

for j=1:size(ContourIndPing_i,1)
    
d=ContourIndPing_i(j,:);
d=d(~isnan(d));
e=ContourIndDepth_i(j,:);
e=e(~isnan(e));
in = inpolygon(TSIndPing70,TSIndDepth70,d,e)
plot(d,e)
hold on
plot(TSIndPing70(in),TSIndDepth70(in),'linestyle','none','marker','+','MarkerEdgeColor',C{i}) % points inside

TSIndPing70(in)
trackLabel70(in)
uniquetrackLabel70 = unique(trackLabel70(in));
TS_mean70 = cell(max(LabelIndex),0);

for k=1:length(uniquetrackLabel70)
    temp = TScomp70(trackLabel70==uniquetrackLabel70(k));
    plot(TSIndPing70(in),TSIndDepth70(in),'linestyle','-','marker','+','MarkerEdgeColor',C{i}) % points inside
    TS_mean70{k} =[TS_mean70 10*log10(nanmean(10.^(temp/10)))];   
end

end
end
set(gca, 'YDir','reverse')




in = inpolygon(TSIndPing70',TSIndDepth70',ContourIndPing,ContourIndDepth)

figure

plot(d,e) % polygon
hold on
plot(TSIndPing70(in),TSIndDepth70(in),,'color',C{ii},'marker','+') % points inside
plot(TSIndPing70(~in),TSIndDepth70(~in),'bo') % points outside
hold off
set(gca, 'YDir','reverse')
