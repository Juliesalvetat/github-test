EI_file='C:\Users\jsalveta\Desktop\Data\FAROFA_2\EK80\HAC\Cruise_FAROFA_2\Treatment20200115_112215\CleanResults\Echointegration\EI001 -Cut\Echointegration.mat';
EI = matfile(EI_file,'Writable',true);
Mask_clean = EI.Mask_clean;
depth_surface = EI.depth_surface;
ind = find(depth_surface>6,1,'first')
Mask_clean(1:ind,:,:) = NaN;
EI.Mask_clean = Mask_clean;

pathwrite=[UserParam.DirSave,'Echointegration.mat']
SaSumCalculation(ei,pathwrite);

ei.Mask_clean(1:8,:, :) = NaN
pathwrite=[UserParam.DirSave,'Echointegration.mat']
SaSumCalculation(ei,pathwrite);