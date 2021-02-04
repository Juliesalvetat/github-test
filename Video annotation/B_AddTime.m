% FAROFA 1
clear
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\Planilhas ajeitadas\FAROFA_1\CSV\';
files = dir([PathName,'**\Processed_G*.csv']);
% Extract only those that are directories.
folders = files.folder;
load('.\data for AddTime\DateTime_Gopro_FAROFA_1.mat','DateTimeStartGopro')
% ajout ST_01 etST_02
for i=1:length(files)
    subdir = dir([files(i).folder,'\Processed_G*.csv']);
    for k=1:length(subdir)
        filename = [subdir(k).folder,'\',subdir(k).name];
        filenewname = [subdir(k).folder,'\With_Time_',subdir(k).name];
        if k == 1
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            number = cell2mat(T.VideoNumber(1));
            number = str2num(number(4:5));
            [Y, M, D, H, MN, S] = datevec(DateTimeStartGopro(number,12:16),'HH:MM');
            StartSec = H*3600+MN*60+S;
            T.Time_UTC = T.Time+StartSec;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,[filenewname])
        else
            
            Time=T.Time_UTC(end);
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            T.Time_UTC = T.Time+Time;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,filenewname)
        end
    end
end
% FAROFA 2 Towed Video
clear
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\Planilhas ajeitadas\FAROFA_2\CSV\TowedVideos\';
files = dir([PathName,'**\Processed_G*.csv']);
%load('C:\Users\garyr\Desktop\JulieHD\Cris\DateTime_Gopro_FAROFA_1.mat')
load('.\data for AddTime\DateTime_VP_DV_FAROFA_2.mat','DateTimeStartDV')
DateTimeStartDV(16:20,:)=DateTimeStartDV(15:19,:);
DateTimeStartDV(15,:) = '20/04/2018 16:40';

for i=1:length(files)
    subdir = dir([files(i).folder,'\Processed_G*.csv']);
    for k=1:length(subdir)
        filename = [subdir(k).folder,'\',subdir(k).name];
        filenewname = [subdir(k).folder,'\With_Time_',subdir(k).name(1:end-4),'.csv'];
        if k ==1
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            number = cell2mat(T.VideoNumber(1));
            number = str2num(number(4:5));
            [Y, M, D, H, MN, S] = datevec(DateTimeStartDV(number,12:16),'HH:MM');
            StartSec = H*3600+MN*60+S;
            T.Time_UTC = T.Time+StartSec;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,[filenewname])
        else
            Time=T.Time_UTC(end);
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            T.Time_UTC = T.Time+Time;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,filenewname)
        end
    end
end


% FAROFA 2 Vertical Profiles
clear
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\FAROFA_2\CSV\VerticalProfiles\';
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\Planilhas ajeitadas\FAROFA_2\CSV\SurfaceTransducer\';

files = dir([PathName,'**\Processed_G*.csv']);
%load('C:\Users\garyr\Desktop\JulieHD\Cris\DateTime_Gopro_FAROFA_1.mat')
load('.\data for AddTime\DateTime_VP_DV_FAROFA_2.mat','DateTimeStartVP')

for i=1:length(files)
    subdir = dir([files(i).folder,'\Processed_G*.csv']);
    for k=1:length(subdir)
        filename = [subdir(k).folder,'\',subdir(k).name];
        filenewname = [subdir(k).folder,'\With_Time_',subdir(k).name(1:end-4),'.csv'];
        if k ==1
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            number = cell2mat(T.VideoNumber(1));
            number = str2num(number(4:5));
            [Y, M, D, H, MN, S] = datevec(DateTimeStartVP(number,12:16),'HH:MM');
            StartSec = H*3600+MN*60+S;
            T.Time_UTC = T.Time+StartSec;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            
            writetable(T,[filenewname])
        else
            
            Time=T.Time_UTC(end);
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            T.Time_UTC = T.Time+Time;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,filenewname)
        end
    end
end
% load('G:\ordi Gary\Area de trabalho\JulieHD/Code/fish/FishExtraction/FAROFA1.mat')

%FAROFA3 Towed Video
clear
%import data importRecapFarofa3.m
path = '.\data for AddTime\RecapVid_FAROFA_3short.csv'; 
opts = detectImportOptions(path);
LogBook = readtable(path,opts);
% select only video instrument from logBook
ix = ismember(LogBook.Instrument,'DV');
Table = table2array(LogBook(ix,4));
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\Planilhas ajeitadas\FAROFA_3\CSV\TowedVideos\';
files = dir([PathName,'**\Processed_G*.csv']);




for i=72:73%length(files)
    subdir = dir([files(i).folder,'\Processed_G*.csv']);
    for k=1:length(subdir)
        filename = [subdir(k).folder,'\',subdir(k).name];
        filenewname = [subdir(k).folder,'\With_Time_',subdir(k).name(1:end-4),'.csv'];
        if k ==1
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            number = cell2mat(T.VideoNumber(1));
            number = str2num(number(4:5));
            [Y, M, D, H, MN, S] = datevec(Table(number),'HH:MM:SS');
            StartSec = H*3600+MN*60+S;
            T.Time_UTC = T.Time+StartSec;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            
            writetable(T,[filenewname])
        else
            
            Time=T.Time_UTC(end);
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            T.Time_UTC = T.Time+Time;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,filenewname)
        end
    end
end


% FAROFA3 ROV
clear
%import data importRecapFarofa3.m
path = '.\data for AddTime\RecapVid_FAROFA_3short.csv'; %FAROFA2
opts = detectImportOptions(path);%FAROFA2
LogBook = readtable(path,opts);


% select only video instrument from logBook
ix = ismember(LogBook.Instrument,'RV');
Table = table2array(LogBook(ix,4));
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\FAROFA_3\CSV\RovProfiles\';
files = dir([PathName,'**\Processed_*.csv']);

for i=1:length(files)
    subdir = dir([files(i).folder,'\Processed_*.csv']);
    for k=1:length(subdir)
        filename = [subdir(k).folder,'\',subdir(k).name];
        filenewname = [subdir(k).folder,'\With_Time_',subdir(k).name(1:end-4),'.csv'];
        if k ==1
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            number = cell2mat(T.VideoNumber(1));
            number = str2num(number(4:5));
            [Y, M, D, H, MN, S] = datevec(Table(number),'HH:MM:SS');
            StartSec = H*3600+MN*60+S;
            T.Time_UTC = T.Time+StartSec;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            
            writetable(T,[filenewname])
        else
            
            Time=T.Time_UTC(end);
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            T.Time_UTC = T.Time+Time;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,filenewname)
        end
    end
end



% FAROFA3 Stereo video
clear
%import data importRecapFarofa3.m
path = '.\data for AddTime\RecapVid_FAROFA_3short.csv'; %FAROFA2
opts = detectImportOptions(path);%FAROFA2
LogBook = readtable(path,opts);


% select only video instrument from logBook
ix = ismember(LogBook.Instrument,'SV');
Table = table2array(LogBook(ix,4));
Table = Table(2:12,:);% on enleve SV_00
PathName = 'C:\Users\jsalveta\Desktop\Planilhas ajeitadas - Copie\Planilhas ajeitadas\FAROFA_3\CSV\StereoVideo\';
files = dir([PathName,'**\Processed_*.csv']);

for i=2:length(files)
    subdir = dir([files(i).folder,'\Processed_*.csv']);
    for k=1:length(subdir)
        filename = [subdir(k).folder,'\',subdir(k).name];
        filenewname = [subdir(k).folder,'\With_Time_',subdir(k).name(1:end-4),'.csv'];
        if k ==1
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            number = cell2mat(T.VideoNumber(1));
            number = str2num(number(4:5));
            [Y, M, D, H, MN, S] = datevec(Table(number),'HH:MM:SS');
            StartSec = H*3600+MN*60+S;
            T.Time_UTC = T.Time+StartSec;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            
            writetable(T,[filenewname])
        else
            
            Time=T.Time_UTC(end);
            opts = detectImportOptions(filename);
            opts = setvartype(opts,'Time',{'double'});
            T = readtable(filename,opts);
            T.Time_UTC = T.Time+Time;
            T.Date=datestr(T.Date, 'yyyy/mm/dd');
            writetable(T,filenewname)
        end
    end
end
