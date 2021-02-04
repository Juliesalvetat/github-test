clear 
%path = 'C:\Users\garyr\Desktop\logbookFAROFA3-2019-05-02.xlsx'; %FAROFA3
%opts = detectImportOptions(path, 'NumHeaderLines', 7); %FAROFA3

%path = 'F:\FAROFA_2\logbookFAROFA2-2019-05-03.xlsx'; %FAROFA2
path = 'C:\Users\jsalveta\Desktop\processing video\RecapVid\logbookFAROFA3-2019-05-02.xlsx'; %FAROFA2


opts = detectImportOptions(path, 'NumHeaderLines', 7);%FAROFA2
opts.VariableNames
opts.VariableTypes 
opts = setvartype(opts,'Time_UTC_',{'datetime'});
LogBook = readtable(path,opts);
LogBook.Date = datestr(LogBook.Date, 'dd/mm/yyyy');
LogBook.Time_UTC_ = datestr(LogBook.Time_UTC_, 'HH:MM');

VideoType = {'SV';'DV';'VP';'RV'};
LogBookShort = [];

% select only video instrument from logBook
for instru=1:length(VideoType)
ix = ismember(LogBook.Instrument,VideoType(instru));
[C,ia,ic] = unique(LogBook.Instrument_number(ix,:));
L = LogBook.Instrument_number(ix,:);
if ~isempty(C)
    LogBookVideoType=LogBook(ix,{'Instrument_number','Instrument','Date','Time_UTC_'});
LogBookShort=[LogBookShort; LogBookVideoType(ia,:)];
end
end

RecapVid=readtable('F:\Cris\RecapVid\RecapVid_FAROFA_2.csv');
VideoType=unique(RecapVid.VideoType);
RecapVidShort=[];
for i=1:size(VideoType,1) % for each Video Type we will calculate the duration 
ix=ismember(RecapVid.VideoType,VideoType(i)); 
RecapVidShort = [RecapVidShort ; RecapVid(ix,[4,7,8,10])] ;
end
[C,ia,ic] = unique(RecapVidShort.VideoNumber,'stable');
RecapVidShort = RecapVidShort(ia,:);
FinalTable=[LogBookShort RecapVidShort(:,[2,3,4])]; % FAROFA2&3
%writetable(FinalTable,['RecapRapport_FAROFA_1.csv'])
 writetable(FinalTable,['RecapRapport_FAROFA_2.csv'])
% writetable(FinalTable,['RecapRapport_FAROFA_3.csv'])

