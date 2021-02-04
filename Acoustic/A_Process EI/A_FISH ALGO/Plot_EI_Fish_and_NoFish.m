clear all
addpath('./privat');
%original EI file :

EI_file='C:\Users\jsalveta\Desktop\DV2 FAROFA2\Cruise_CruiseName\Treatment20191122_101152\CleanResults\Echointegration\EI001\Echointegration.mat';
%EI_file='G:\FAROFA1\EK80\HAC\Cruise_ABRACINHOS01\Treatment20171207_101445\SubTreat20180313_215604_-90_1ping_02\CleanResults\Echointegration\Echointegration.mat';
%EI_file='G:\Processing_2019\EI_cut_fish\3zones.mat';

EI = matfile(EI_file,'Writable',false);
frequ_to_plot=2;

%Echointegration file containing only classified data
%Class_Filename='C:\Gildas\Abrazos_2\Processing_2019\Fish_Algo\Echointegration_FDN_ADR_FISH';
%Class_Filename='G:\FAROFA1\EK80\HAC\Cruise_ABRACINHOS01\Treatment20171207_101445\SubTreat20180313_215604_-90_1ping_02\CleanResults\Echointegration\Echointegration_fish_no_fish.mat';
%Class_Filename='G:\Processing_2019\EI_cut_fish\3zones_fish_no_fish.mat';
Class_Filename='C:\Users\jsalveta\Desktop\DV2 FAROFA2\Cruise_CruiseName\Treatment20191122_101152\CleanResults\Echointegration\EI001\Echointegration_FISH\Echointegration.mat';

Classif = matfile(Class_Filename,'Writable',false);
depth_vect=Classif.depth_surface;
Bottom=Classif.Bottom;
maskFish=Classif.Mask_clean;
SvFish=Classif.Sv_surface;
SvFish=SvFish.*maskFish;

SvNoFish=EI.Sv_surface;
maskNofish=ones(size(maskFish));
maskNofish(maskFish==1)=NaN;
SvNoFish=SvNoFish.*maskNofish;
Sv_surface=Classif.Sv_surface;
figure()
Ping_to_plot=[1:1200];
subplot(2,1,2)
gtsecho(SvFish(:,Ping_to_plot,2),depth_vect,[-80 -40],Bottom(1,Ping_to_plot));
ylim([0 max(Bottom(1,Ping_to_plot))])
 title('Fish algorithm result - FAROFA 2 DV2 ')
subplot(2,1,1)
%figure
gtsecho(Sv_surface(:,Ping_to_plot,2),depth_vect,[-80 -40],Bottom(1,Ping_to_plot));
ylim([0 max(Bottom(1,Ping_to_plot))])
 title('Original Echointegration- FAROFA 2 DV2 ')

ylim([0 110])
% title('No Fish algorithm result - FDN ')
title('Original 70kHz')
% gtsecho(SvNoFish(:,Ping_to_plot,frequ_to_plot),depth_vect,[-80 -40],Bottom(frequ_to_plot,Ping_to_plot));
% ylim([0 110])
% title('No Fish algorithm result - FDN ')