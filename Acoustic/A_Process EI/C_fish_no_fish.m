addpath("C:\Users\jsalveta\Desktop\Matlab code\gfunctions")
% load EI FISH
filename = '.\Echointegration_FISH_cleaned\Echointegration.mat';
EIfish = matfile(filename,'Writable',false); clear filename;
depth_surface = EIfish.depth_surface;
Bottom = EIfish.Bottom(1,:);
Sv_surface_fish_70 = EIfish.Sv_surface(:,:,1);
Sv_surface_fish_200 = EIfish.Sv_surface(:,:,2);
Mask_clean = EIfish.Mask_clean(:,:,1);
Sv_surface_fish_70 = Sv_surface_fish_70.*Mask_clean;
Sv_surface_fish_200 = Sv_surface_fish_200.*Mask_clean;

% figure;gtsecho(Sv_surface_fish_70(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('Fish EI 70 kHz');ylim([0 100])
% figure;gtsecho(Sv_surface_fish_200(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('Fish EI 200 kHz');ylim([0 100])

%load EI Originale
filename = 'Echointegration_original\Echointegration.mat';
EIorigin = matfile(filename,'Writable',false); clear filename;
Sv_surface_origin_70 = EIorigin.Sv_surface(:,:,1);
Sv_surface_origin_200 = EIorigin.Sv_surface(:,:,2);
Sa_surface_origin_70 = EIorigin.Sa_surface(:,:,1);
Sa_surface_origin_200 = EIorigin.Sa_surface(:,:,2);
% applique Mask_clean (correction manuelles NaN)
Mask_clean_origin = EIorigin.Mask_clean(:,:,1);
Sv_surface_origin_70  = Sv_surface_origin_70.*Mask_clean_origin;
Sv_surface_origin_200 = Sv_surface_origin_200.*Mask_clean_origin;
Sa_surface_origin_70  = Sa_surface_origin_70.*Mask_clean_origin;
Sa_surface_origin_200 = Sa_surface_origin_200.*Mask_clean_origin;

%  figure;gtsecho(Sv_surface_origin_70(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('Original EI 70 kHz');ylim([0 100])
% figure;gtsecho(Sv_surface_origin_200(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('Original EI 200 kHz');ylim([0 100])

% Mask no fish
mask_no_fish=isnan(Sv_surface_fish_200);
no_fish_200=gmask(Sv_surface_origin_200,mask_no_fish,-150);
no_fish_70=gmask(Sv_surface_origin_70,mask_no_fish,-150);
sa_no_fish_200=gmask(Sa_surface_origin_200,mask_no_fish,0);
sa_no_fish_70=gmask(Sa_surface_origin_70,mask_no_fish,0);

% figure;gtsecho(no_fish_70(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('Sv no fish 70 kHz');ylim([0 100])
% figure;gtsecho(no_fish_200(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('Sv no fish 200 kHz');ylim([0 100])
% figure;gtsecho(sa_no_fish_70(:,57312:59150),depth_surface',[0 40],Bottom(57312:59150));colorbar;title('Sa fish 70 kHz');ylim([0 100])
% figure;gtsecho(sa_no_fish_200(:,57312:59150),depth_surface',[0 40],Bottom(57312:59150));colorbar;title('Sa no fish EI 200 kHz');ylim([0 100])

%recupere des leftover fish
mask_leftover_fish = no_fish_70>-55 & no_fish_200>-55;%| (no_fish_70>-60 & mask_75m) & no_fish_200>-60;
leftover_fish_200=gmask(no_fish_200,mask_leftover_fish,-150);
leftover_fish_70=gmask(no_fish_70,mask_leftover_fish,-150);
figure;gtsecho(leftover_fish_70(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('leftover fish 70 kHz');ylim([0 100])
figure;gtsecho(leftover_fish_200(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('leftover fish 200 kHz');ylim([0 100])

% offset 
maskBottom = zeros(size(no_fish_70));
for i=1:length(Bottom)
    [minValue,closestIndex2] = min(abs(depth_surface-(Bottom(i)-0.75)));
    maskBottom(closestIndex2+1:end,i)=ones(size(maskBottom,1)-closestIndex2,1);
end
figure;imagesc(maskBottom(:,57312:59150));colorbar
maskoffset = maskBottom & (no_fish_70>-55 | no_fish_200>-55);
figure;imagesc(maskoffset(:,57312:59150));colorbar
figure;imagesc(maskoffset(:,6947:7214));colorbar

mask_fish_final = mask_leftover_fish  | ~mask_no_fish | maskoffset;
figure;imagesc(mask_fish_final(:,6947:7214));colorbar
figure;imagesc(mask_fish_final(:,57312:59150));colorbar

fish_200_final=gmask(Sv_surface_origin_200,mask_fish_final,-150);
fish_70_final=gmask(Sv_surface_origin_70,mask_fish_final,-150);
sa_fish_200_final=gmask(Sa_surface_origin_200,mask_fish_final,0);
sa_fish_70_final=gmask(Sa_surface_origin_70,mask_fish_final,0);

figure;gtsecho(fish_70_final(:,6947:7214),depth_surface',[-90 -40],Bottom(6947:7214));colorbar;title('fish final 70 kHz');ylim([0 100])
figure;gtsecho(fish_70_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 70 kHz');ylim([0 100])
figure;gtsecho(fish_200_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 200 kHz');ylim([0 100])


filename = '.\FISH and leftover FISH\Echointegration.mat';
EIfish = matfile(filename,'Writable',true); clear filename;
EIfish.Sv_surface(:,:,1)=fish_70_final;
EIfish.Sv_surface(:,:,2)=fish_200_final;
EIfish.Sa_surface(:,:,1)=sa_fish_70_final;
EIfish.Sa_surface(:,:,2)=sa_fish_200_final;


filename = '.\FISH -Manual corrections\Echointegration.mat';
EIfish = matfile(filename,'Writable',true); clear filename;
depth_surface = EIfish.depth_surface;
Bottom = EIfish.Bottom(1,:);
Mask_clean = EIfish.Mask_clean(:,:,1);
Mask = isnan(Mask_clean);
%figure;imagesc(Mask(:,57312:59150))-
Sv_surface_70 = EIfish.Sv_surface(:,:,1);
%figure;gtsecho(Sv_surface_70(:,57312:59150),depth_surface',[90 -40],Bottom(57312:59150));colorbar;title('fish final 70 kHz');ylim([0 100])

Mask = isnan(Mask_clean) | isnan(Sv_surface_70) ;
Mask = ~ Mask;
%figure;imagesc(Mask(:,57312:59150))

Sv_surface_200 = EIfish.Sv_surface(:,:,2);
Sa_surface_70 = EIfish.Sa_surface(:,:,1);
Sa_surface_200 = EIfish.Sa_surface(:,:,2);

fish_200_final=gmask(Sv_surface_200,Mask,-150);
fish_70_final=gmask(Sv_surface_70,Mask,-150);
sa_fish_200_final=gmask(Sa_surface_200,Mask,0);
sa_fish_70_final=gmask(Sa_surface_70,Mask,0);

figure;gtsecho(fish_70_final(:,826184:831081),depth_surface',[-90 -40],Bottom(826184:831081));colorbar;title('fish final 70 kHz');ylim([0 100])
% figure;gtsecho(fish_200_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 200 kHz');ylim([0 100])

filename = '.\Echointegration_original\Echointegration.mat';
EIorigin = matfile(filename,'Writable',false); clear filename;
Mask_clean_origin = EIorigin.Mask_clean(:,:,1);
fish_70_final  = fish_70_final.*Mask_clean_origin;
fish_200_final = fish_200_final.*Mask_clean_origin;
sa_fish_70_final  = sa_fish_70_final.*Mask_clean_origin;
sa_fish_200_final = sa_fish_200_final.*Mask_clean_origin;

% figure;gtsecho(fish_70_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 70 kHz');ylim([0 100])
% figure;gtsecho(fish_200_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 200 kHz');ylim([0 100])

EIfish.Sv_surface(:,:,1)=fish_70_final;
EIfish.Sv_surface(:,:,2)=fish_200_final;
EIfish.Sa_surface(:,:,1)=sa_fish_70_final;
EIfish.Sa_surface(:,:,2)=sa_fish_200_final;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no fish
filename = '.\FISH -Manual corrections\Echointegration.mat';
EIfish = matfile(filename,'Writable',true); clear filename;
depth_surface = EIfish.depth_surface;
Bottom = EIfish.Bottom(1,:);
Sv_surface_70 = EIfish.Sv_surface(:,:,1);
Mask = Sv_surface_70==-150;
% figure;gtsecho(Sv_surface_70(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 70 kHz');ylim([0 100])
% figure;imagesc(Mask(:,57312:59150))
% figure;imagesc(Mask(:,6947:7214))
% figure;gtsecho(Sv_surface_70(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('fish final 200 kHz');ylim([0 100])
% figure;gtsecho(Sv_surface_70(:,6947:7214),depth_surface',[-90 -40],Bottom(6947:7214));colorbar;title('fish final 70 kHz');ylim([0 100])
% figure;imagesc(~Mask(:,6947:7214))
% figure;imagesc(~Mask(:,57312:59150))


filename = '.\NO FISH\Echointegration.mat';
EInofish = matfile(filename,'Writable',true); clear filename;

no_fish_200_final=gmask(Sv_surface_origin_200,Mask,-150);
no_fish_70_final=gmask(Sv_surface_origin_70,Mask,-150);
sa_no_fish_200_final=gmask(Sa_surface_origin_200,Mask,0);
sa_no_fish_70_final=gmask(Sa_surface_origin_70,Mask,0);

% figure;gtsecho(no_fish_70_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('no fish final 70 kHz');ylim([0 100])
% figure;gtsecho(no_fish_200_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('no fish final 200 kHz');ylim([0 100])
% figure;gtsecho(no_fish_70_final(:,6947:7214),depth_surface',[-90 -40],Bottom(6947:7214));colorbar;title('no fish final 70 kHz');ylim([0 100])
% figure;gtsecho(no_fish_200_final(:,57312:59150),depth_surface',[-90 -40],Bottom(57312:59150));colorbar;title('no fish final 200 kHz');ylim([0 100])
figure;gtsecho(no_fish_70_final(:,826184:831081),depth_surface',[-90 -40],Bottom(826184:831081));colorbar;title('fish final 70 kHz');ylim([0 100])

EInofish.Sv_surface(:,:,1)=no_fish_70_final;
EInofish.Sv_surface(:,:,2)=no_fish_200_final;
EInofish.Sa_surface(:,:,1)=sa_no_fish_70_final;
EInofish.Sa_surface(:,:,2)=sa_no_fish_200_final;
%%%%%%%%%%%%% ORIGINAL integration mask clean
clear all; close all; clc
folder = {'FISH' 'NO FISH' 'ORIGINAL'};

for FAROFA=1:3
            
         for f=2:length(folder)
        
        Class_Filename = ['.\FAROFA',num2str(FAROFA),'\',folder{f},'\Echointegration.mat'];
                Classif = matfile(Class_Filename,'Writable',true);

% Bottom = Classif.Bottom(1,:);
% depth_surface = Classif.depth_surface;
Sv_surface = Classif.Sv_surface;
Sa_surface = Classif.Sa_surface;
Mask_clean = Classif.Mask_clean;
Classif.Sa_surface = Sa_surface.*Mask_clean;
Classif.Sv_surface = Sv_surface.*Mask_clean;
         end
end
%Mask_clean(isnan(Mask_clean))=0;
%Mask_clean = logical(Mask_clean);
% Classif.Sv_surface(:,:,2)=gmask(Sv_surface_fish(:,:,2),Mask_clean(:,:,2),-150);
% Classif.Sa_surface(:,:,1)=gmask(Sv_surface_fish(:,:,1),Mask_clean(:,:,1),-150);

figure;
subplot(4,1,1:3)
gtsecho(Classif.Sa_surface(:,315282:315799,1),depth_surface',[0 40],Bottom(315282:315799));colorbar;title('Fish EI 70 kHz');ylim([0 100])
subplot(4,1,4)
plot(nansum(Classif.Sa_surface(:,315282:315799,1),1))

figure;
subplot(4,1,1:3)
gtsecho(Classif.Sa_surface(:,315282:315799,2),depth_surface',[0 40],Bottom(315282:315799));colorbar;title('Fish EI 70 kHz');ylim([0 100])
subplot(4,1,4)
plot(nansum(Classif.Sa_surface(:,315282:315799,2),1))

%%%%%%%%%% Sa 
clear all; close all; clc

folder = {'FISH' 'NO FISH' 'ORIGINAL'};

for FAROFA=1:3
    
    for f=1:length(folder)
        disp(['FAROFA ',num2str(FAROFA),' : ',folder{f}])

        Class_Filename = ['.\FAROFA',num2str(FAROFA),'\',folder{f},'\Echointegration.mat'];        
        Classif = matfile(Class_Filename,'Writable',true);

Bottom = Classif.Bottom(1,:);
depth_surface = Classif.depth_surface;
Bottom = Classif.Bottom(1,:);
Sv_surface_70 = Classif.Sv_surface(:,:,1);
Sv_surface_200 = Classif.Sv_surface(:,:,2);
Sa_surface_70 = Classif.Sa_surface(:,:,1);
Sa_surface_200 = Classif.Sa_surface(:,:,2);


% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_200(:,315282:315799),depth_surface',[0 40],Bottom(315282:315799));colorbar;title('Fish EI 200 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_200(:,315282:315799),1))
% 
% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_70(:,315282:315799),depth_surface',[0 40],Bottom(315282:315799));colorbar;title('Fish EI 70 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_70(:,315282:315799),1))


maskBottom = zeros(size(Sa_surface_70));
for i=1:length(Bottom)
    [minValue,closestIndex2] = min(abs(depth_surface-(Bottom(i)-1.5)));
    maskBottom(closestIndex2+1:end,i)=ones(size(maskBottom,1)-closestIndex2,1);
    %SaSum_Bottom(i) = nansum(
end
%figure;imagesc(maskBottom(:,315282:315799));colorbar
maskoffset_Sa = ~(maskBottom & (Sa_surface_70>5000 | Sa_surface_200>1000));
%figure;imagesc(maskoffset_Sa(:,315282:315799));colorbar

Sv_surface_200_final=gmask(Sv_surface_200,maskoffset_Sa,NaN);
Sv_surface_70_final=gmask(Sv_surface_70,maskoffset_Sa,NaN);
Sa_surface_200_final=gmask(Sa_surface_200,maskoffset_Sa,NaN);
Sa_surface_70_final=gmask(Sa_surface_70,maskoffset_Sa,NaN);

% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_200_final(:,315282:315799),depth_surface',[0 40],Bottom(315282:315799));title('Fish EI 200 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_200_final(:,315282:315799),1)); xlim([0 length(315282:315799)+1])
% 
% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_70_final(:,315282:315799),depth_surface',[0 40],Bottom(315282:315799));title('Fish EI 70 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_70_final(:,315282:315799),1));xlim([0 length(315282:315799)+1])
% 
% 
% 
% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_200_final(:,551805:552805),depth_surface',[0 40],Bottom(551805:552805));title('Fish EI 200 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_200_final(:,551805:552805),1)); xlim([0 length(551805:552805)+1])
% 
% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_70_final(:,551805:552805),depth_surface',[0 40],Bottom(551805:552805));title('Fish EI 200 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_70_final(:,551805:552805),1)); xlim([0 length(551805:552805)+1])
% 
% 
% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_200(:,551805:552805),depth_surface',[0 40],Bottom(551805:552805));title('Fish EI 200 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_200(:,551805:552805),1)); xlim([0 length(551805:552805)+1])
% 
% figure;
% subplot(4,1,1:3)
% gtsecho(Sa_surface_70_final(:,551805:552805),depth_surface',[0 40],Bottom(551805:552805));title('Fish EI 200 kHz');ylim([0 100])
% subplot(4,1,4)
% plot(nansum(Sa_surface_70_final(:,551805:552805),1)); xlim([0 length(551805:552805)+1])

Classif.Sv_surface(:,:,1) = Sv_surface_70_final;
Classif.Sv_surface(:,:,2) = Sv_surface_200_final;

Classif.Sa_surface(:,:,1) = Sa_surface_70_final;
Classif.Sa_surface(:,:,2) = Sa_surface_200_final;

    end    
end