fclose all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read user parameters for processing
FISH_ALGO_Define_param_forMatecho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./privat')
%addpath('..\Zscore\codes_diff_and_sum\Tool_HacGenericProcess\nansuite');
%dataPath = 'G:\FAROFA2\EK80\HAC_transects_FdN\Cruise_transectsFdN\Treatment20180803_101619\SubTreat20190111_135138_-90_1ping_02\CleanResults\Echointegration\Echointegration.mat';
%dataPath = 'G:\FAROFA1\EK80\HAC\Cruise_ABRACINHOS01\Treatment20171207_101445\SubTreat20180313_215604_-90_1ping_02\CleanResults\Echointegration\Echointegration.mat';
%dataPath = 'G:\FAROFA2\EK80\HAC_peche\Cruise_peche23avril_Julie\Treatment20180824_145132\SubTreat20190111_091424_-90_1ping_02\CleanResults\Echointegration\Echointegration.mat';
%dataPath = 'G:\FAROFA2\EK80\HAC_transects_AFD\Cruise_transectsAFD_Anne\Treatment20180426_135008\CleanResults\Echointegration\Echointegration.mat';
%outputFilename = 'G:\Processing_j\Echointegration\FAROFA2\seuilSv-65_seuilSum-120\peche\Echointegration.mat';

%copy data in new file name
input_echo = matfile(dataPath,'Writable',false);
echo = matfile(outputFilename,'Writable',true);

copyfile(dataPath,outputFilename)

FrequencySort = input_echo.FrequencySort;
idF01=find(FrequencySort/1000==Freq(1));
idF02=find(FrequencySort/1000==Freq(2));
idF03=find(FrequencySort/1000==Freq(3));

depth_surface=echo.depth_surface';
id_maxD=find(depth_surface<Max_depth,1,'last');
nlayer_Toadd=length(depth_surface)-id_maxD;
Night1Sunrise2Day3Sunset4=echo.Night1Sunrise2Day3Sunset4;


%we work on all ESU...
ESU_SELection=[1,size(input_echo.Time,2)];
Loop_sel_esu=ESU_SELection;

if ESU_SELection(end)-ESU_SELection(1) > ESU_max_N
    
    sequence_deb=ESU_SELection(1):ESU_max_N:ESU_SELection(end);
    if  sequence_deb(end)+Size_filt_horiz>ESU_SELection(end)
    sequence_deb = sequence_deb(1:end-1);
    end
    Loop_sel_esu=NaN(length(sequence_deb),2);
    Loop_sel_esu(1,1)=ESU_SELection(1);
    Loop_sel_esu(1,2)=sequence_deb(2)+Size_filt_horiz;
    Loop_sel_esu(2:length(sequence_deb),1)=sequence_deb(2:end)-Size_filt_horiz;
    Loop_sel_esu(2:length(sequence_deb)-1,2)=sequence_deb(2:end-1)+ESU_max_N+Size_filt_horiz;
    Loop_sel_esu(end,2)=ESU_SELection(end);
end

nloop=size(Loop_sel_esu,1);

h = waitbar(0,'Please wait...');


for i_loop=1:nloop
    disp(['Loop #',num2str(i_loop),' over ', num2str(nloop)])
    
    %waitbar(i_loop/nloop,h)
    
    %update ESU selection
    ESU_SELection=[Loop_sel_esu(i_loop,1):Loop_sel_esu(i_loop,2)];
    
    
    Loop_algo
    
    %remove esu at beginning and end of selection (that was selected
    %according to horizontal size of filter)
    if i_loop==1 %if first loop
        if size(mask_fish,1)<=ESU_max_N %if only one loop and number of esu is less than N esu in loop
            size_to_use=1:size(mask_fish,1);
            concat_fish
        else
            size_to_use=1:size(mask_fish,1)-(Size_filt_horiz+1);
            concat_fish
        end
    elseif i_loop==nloop %if last loop
        size_to_use=Size_filt_horiz+1:size(mask_fish,1);
        concat_fish
        
    else %if not first or last loop
        size_to_use=Size_filt_horiz+1:size(mask_fish,1)-(Size_filt_horiz+1);
        concat_fish
    end
    
    if Process_stop_loop %if pause is required ...ask an input from user
        input('Type Enter to continue\n','s');
    end
    
end %end Loop by esus sequences

%Add zeros to Cat_Fish_mask deaper than Max_depth used for algorithm



clear Cat_Fish_mask
%No_Fish_mask=~mask_fish;

%ESU selection update before saving
%ESU_SELection=[Loop_sel_esu(1,1):Loop_sel_esu(end,2)];

%Save masks


%clear temp
%save(outputFilename,'No_Fish_mask','mask_fish','dataPath','Max_depth','ESU_SELection','outputFilename')

%delete .asv files
delete('*.asv');
delete('./privat/*.asv');

%uncomment to check whole ESU selection
%    Loop_algo











