%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters to update by user for step D (echogram for three frequencies RGB color scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path to EI from matecho)
%file_path ='C:\Users\garyr\Desktop\JulieHD\cruze\Test_TS\Cruise_CruiseName_newTest\Treatment20190509_175138\CleanResults\Echointegration\';
[baseFileName, folder] = uigetfile(pwd,'Select Echointegration.mat');
dataPath=[folder,baseFileName];


% Create Echointegration_FISH folder if it doesn't already exist
if ~exist([folder, 'Echointegration_FISH'], 'dir')
  mkdir([folder, '\Echointegration_FISH']);
end

outputFilename=[folder, '\Echointegration_FISH\Echointegration.mat'];

%Maximum number of ESU in each processing loop.
ESU_max_N=1000;

%Processing stop at each loop  ? set to one to say YES
Process_stop_loop=0;

%The Three Frequencies (in kHz) used to plot echogram with RGB color scale
%First one will be Red, second one Green, and third one Blue
Freq=[70 200 200]; %in kHz

%max depth in m to use in echogram
Max_depth=110; %in m

%low limits (in dB) for each frequency. All Sv
%values under this limits will be convert to zero. If all frequecies are set to zero, pixel will be black.
%User will be able to adjust this values by clicking on figure when plot
%will be done.
Low_Limits=[-80 -80 -80]; %in dB

%Range (in dB) from low limits for each frequency. Sv values greater than
%Low_Limits+Ranges will be set to 255. If all frequencies are set to 255, pixel will be white.
Ranges=[40 40 40]; %in dB

%option to do convolution on echogram
Do_convolution=0;
%Size of square convolution filter (size in number of horizontal and vertical samples) 
Size_conv=1;


%%%%%%%%%%%%%%%%%%%%%%% ALGO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Half Filter size to compute variance
%If Size_filt_horiz=10, variance will be compute 10 esus before and 10 esus
%after current sample. Idem in vertical direction
Size_filt_vertical=2; %in number of vertical samples
Size_filt_horiz=15; %in number of esu

%If frequency 1 and 2 are both below this treshold, EI cell will not be use
Sv_LowThresh=-80;
%Trheshold on Variance for first and second Sv
thresh_var_1=70;
thresh_var_2=90;
%Threshold on Sv for first High variance classification
Tresh_sv1=-65;% in dB
Tresh_sv2=-65;% in dB
%Threshold on Sv sum for very dense fish classification without high
Thresh_sum=-120;


