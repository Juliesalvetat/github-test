%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construsting the table RecapVid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note to my self try doing it in R next time 
% Create a .csv with all video detailed using the arborescence to identify each Video and ffmpeg to get video Metadata such as time and duration 

clear
addpath('E:\ordi Gary\Area de trabalho\JulieHD\Cris\privat\ffmpeg') % Add path to ffmpeg function to collect time and video duration 

%%%  Enter the name of your cruise 
CruiseName = 'FAROFA'; 
%%%  Enter how many cruise you want to create a RecapVid for 
CruiseNumber = 3;

for n=1:CruiseNumber % nomber of cruise
    % define path for each Cruise
    if n== 1 | n== 2 
       PathName=['E:\VIDEO_HD\',CruiseName,'_',num2str(n),'\VIDEOS\']; % FAROFA 1 & 2 are on the same disk 
    elseif n==3
       PathName=['E:\VIDEO_HD\',CruiseName,'_',num2str(n),'\VIDEOS\'];
    end
    % List all .MP4 in folder in PathName and subfolders
    files = dir([PathName,'**\*.MP4']);    
    % Create empty table containing 
    RecapVid  = cell2table(cell(size(files,1),10), 'VariableNames', {'Cruise','VideoType','Date','VideoNumber','VideoName','Processed','ffmpegDate','ffmpegTime','ffmpegDuration_seconds','ffmpegDuration_hhmmss'});
    for i=1:size(files,1)
        C = strsplit(files(i).folder,'\'); % separer utilisant T comme delimitateur
        Cruise (i,:) = C(3);
        VideoType(i,:) = C(5);
        Date(i,:) = C(6);
        VideoNumber(i,:) = C(7);
        VideoName(i,:) = files(i).name;
        % Check if there is csv in folder and subfolders
        if ~isempty(ls([files(i).folder,'\','**\Processed*.csv']))
            Processed(i,:)='yes';
        else
            Processed(i,:)='no';
        end
        INFO = ffmpeginfo([files(i).folder,'\',files(i).name]);
%         ffmpegTime(i,:)=[INFO.meta.creation_time(1:10),' ',INFO.meta.creation_time(12:22)];
        ffmpegDate(i,:) = INFO.meta.creation_time(1:10);
        ffmpegTime(i,:) = INFO.meta.creation_time(12:22);

        % C = strsplit(INFO.meta.creation_time,'T') % other way to do it using T as separator 
        % regexp(INFO.meta.creation_time, '[TZ]', 'split') % other way to do it using T and Z as separator 
        ffmpegDuration_seconds(i,:)=INFO.duration;  % in seconds
        ffmpegDuration_hhmmss(i,:) = datestr(seconds(INFO.duration),'HH:MM:SS');

    end

RecapVid.Cruise = Cruise;
RecapVid.VideoType = VideoType;
RecapVid.Date = Date;
RecapVid.VideoNumber = VideoNumber;
RecapVid.VideoName = VideoName;
RecapVid.Processed = Processed;
RecapVid.ffmpegDate = ffmpegDate;
RecapVid.ffmpegTime = ffmpegTime;
RecapVid.ffmpegDuration_seconds = ffmpegDuration_seconds;
RecapVid.ffmpegDuration_hhmmss = ffmpegDuration_hhmmss;

writetable(RecapVid,['RecapVid_', CruiseName,'_', num2str(n),'.csv'])
clear RecapVid Cruise VideoType Date VideoNumber VideoName Processed ffmpegDate ffmpegTime ffmpegDuration_hhmmss ffmpegDuration_seconds INFO
end

%%%%%%%%%%%%%%%%%%%%%% Read RecapVid and culculate Time duration for each VideoType 
clear 
CruiseName = 'FAROFA'; 
Path2csv = 'F:\Cris\RecapVid\';
for n=1:3
csvName = ['RecapVid_FAROFA_',num2str(n),'.csv'];
RecapVid=readtable([Path2csv,csvName]);
VideoType=unique(RecapVid.VideoType);
for i=1:size(VideoType,1) % for each Video Type we will calculate the duration 
ix=ismember(RecapVid.VideoType,VideoType(i));  % find in column RecapVid.VideoType the Videotype for which we you sum duration
sec=sum(RecapVid.ffmpegDuration_seconds(ix,:)) ;
Number(i,:) = length(unique(RecapVid.VideoNumber(ix,:)));
if isequal(cell2mat(VideoType(i)),'StereoVideo')
Duration(i,:) = datestr(seconds(sec)/2,'HH:MM:SS');
else
    Duration(i,:) = datestr(seconds(sec),'HH:MM:SS');
end
end
RecapVidDuration.Cruise = RecapVid.Cruise(1:size(VideoType,1));
RecapVidDuration.VideoType = VideoType;
RecapVidDuration.TotalNumber = Number;
RecapVidDuration.TotalDuration = Duration;
RecapVidDuration = struct2table(RecapVidDuration);
writetable(RecapVidDuration,['RecapVidDuration_', CruiseName,'_', num2str(n),'.csv'])
clear RecapVid RecapVidDuration Duration VideoType Cruise Number
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Concatenate RecapVid
clear 
CruiseName = 'FAROFA'; 
RecapVid  = [];
RecapVidDuration = [];
Path2csv = 'F:\Cris\';
for n=1:3
csvNameR = ['RecapVid_FAROFA_',num2str(n),'.csv'];
csvNameRD = ['RecapVidDuration_FAROFA_',num2str(n),'.csv'];
R=readtable([Path2csv,csvNameR]); 
RD=readtable([Path2csv,csvNameRD]);
RecapVid=[RecapVid; R];
RecapVidDuration=[RecapVidDuration; RD];
end 
writetable(RecapVid,['RecapVid_all_', CruiseName,'.csv'])
writetable(RecapVidDuration,['RecapVidDuration_all_', CruiseName,'.csv'])

%%%%%