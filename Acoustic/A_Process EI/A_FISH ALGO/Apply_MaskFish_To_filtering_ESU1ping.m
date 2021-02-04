clear all; close all
%clock

%Echointegration file containing only classified data
%Echointegration_FISH_Path='C:\Users\garyr\Desktop\JulieHD\cruze\Test_TS\Cruise_CruiseName_newTest\Treatment20190509_175138\CleanResults\Echointegration\Echointegration_FISH\Echointegration.mat';
[baseFileName, folder] = uigetfile(pwd,'Select Echointegration_FISH_cleaned\Echointegration.mat');
Echointegration_FISH_Path=[folder,baseFileName];
%echogram file containing data of the cruise
%Echogram_Path='C:\Users\garyr\Desktop\JulieHD\cruze\Test_TS\Cruise_CruiseName_newTest\MatData\Echogram.mat';
[baseFileName, folder] = uigetfile(pwd,'Select \MatData\Echogram.mat');
Echogram_Path=[folder,baseFileName];
%filtering.mat of coresponding matecho cruise to update with classified
%data
%Filtering_Path='C:\Users\garyr\Desktop\JulieHD\cruze\Test_TS\Cruise_CruiseName_newTest\Treatment20190510_080133_FISH\DataFilterResults\Filtering.mat';
[baseFileName, folder] = uigetfile(pwd,'Select Treatment_FISH\DataFilterResults\Filtering.mat');
Filtering_Path =[folder,baseFileName];
Classif = matfile(Echointegration_FISH_Path,'Writable',true);
echo = matfile(Echogram_Path,'Writable',true);
filt = matfile(Filtering_Path,'Writable',true);

%save old mask if not alredy done!!!!
names_F=fieldnames(filt);


if ~sum(strcmp(names_F,'Mask_Cleaning70_OLD'))% compare names_F to 'Mask_Cleaning70_OLD' if no 'Mask_Cleaning70_OLD' -> 0
    Esu_step = 1000;
    totalNesu=size(filt.CleanBottom,2);
    seq_ESU=1:Esu_step:totalNesu;
    seq_ESU=[seq_ESU,totalNesu];
    for i= 1:length(seq_ESU)-1
        i
        if i==1 % to first create  filt.Mask_Cleaning200_OLD and  filt.Mask_Cleaning70_OLD
            filt.Mask_Cleaning200_OLD=filt.Mask_Cleaning200(:,seq_ESU(i):seq_ESU(i+1));
            filt.Mask_Cleaning70_OLD=filt.Mask_Cleaning70(:,seq_ESU(i):seq_ESU(i+1));
        else
            filt.Mask_Cleaning200_OLD(:,seq_ESU(i)+1:seq_ESU(i+1))=filt.Mask_Cleaning200(:,seq_ESU(i)+1:seq_ESU(i+1));
            filt.Mask_Cleaning70_OLD(:,seq_ESU(i)+1:seq_ESU(i+1))=filt.Mask_Cleaning70(:,seq_ESU(i)+1:seq_ESU(i+1));
        end
    end
    clear Esu_step seq_ESU totalNesu
    
end

clear names_F

%filt.Mask_Cleaning200=filt.Mask_Cleaning200_OLD; % OUT MEMORIE
%filt.Mask_Cleaning70=filt.Mask_Cleaning70_OLD;

%adjust number of esus to process in each loop. if too big, matlab may
%crash
Esu_step=1000;

%Depth limit of algorithm, deaper than this limit the program will set data
%to NaN
Depth_limit_MaskFish=100; %in m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Esus sequence to process in loops
totalNesu=size(Classif.Bottom,2);
seq_ESU=[1:Esu_step:totalNesu,totalNesu];

depth_EI=Classif.depth_surface';
vert_size_layer=depth_EI(2)-depth_EI(1); %layer thikness of Ei

depth_samples_1=echo.Depth70(:,1);
depth_samples_2=echo.Depth200(:,1);

%astuce pour plus tard...
%Set a samples in echogram to -20db just above first EI layer
ind_sample_20dB_1=find ((depth_samples_1<depth_EI(1)-(vert_size_layer) ) & (depth_samples_1>depth_EI(1)-2*(vert_size_layer) ));
ind_sample_20dB_2=find ((depth_samples_2<depth_EI(1)-(vert_size_layer) ) & (depth_samples_2>depth_EI(1)-2*(vert_size_layer) ));

%we don't need any more echogram.mat.
%clear echo

indEI_forsampl_1=NaN(size(depth_samples_1));
indEI_forsampl_2=NaN(size(depth_samples_2));

EsuPingStartEnd=Classif.EsuPingStartEnd;

for id=1:size(depth_EI,2) %for each EI layer..
    matching_depth_1 = (depth_samples_1 >= (depth_EI(id)-(vert_size_layer/2))) & (depth_samples_1 <= (depth_EI(id)+(vert_size_layer/2)));
    matching_depth_2 = (depth_samples_2 >= (depth_EI(id)-(vert_size_layer/2))) & (depth_samples_2 <= (depth_EI(id)+(vert_size_layer/2)));
    %for each vertical sample we have the corresponding index of depth
    %layer in Ei.
    indEI_forsampl_1(matching_depth_1)=id;
    indEI_forsampl_2(matching_depth_2)=id;
    
end

%last indice of sample used for Ei. Deaper than this, Ei was not done
ind_maxdepth_EI=find(depth_EI<=Depth_limit_MaskFish+vert_size_layer/2,1,'last');

ind_minsampledepth_1=find(depth_samples_1<=depth_EI(1)-vert_size_layer/2,1,'last');
ind_minsampledepth_2=find(depth_samples_2<=depth_EI(1)-vert_size_layer/2,1,'last');

ind_maxsampledepth_1=find(depth_samples_1<=Depth_limit_MaskFish+vert_size_layer/2,1,'last');
ind_maxsampledepth_2=find(depth_samples_2<=Depth_limit_MaskFish+vert_size_layer/2,1,'last');%loop by Esu sequences

n_depth_EI=find(depth_EI<=Depth_limit_MaskFish+vert_size_layer/2,1,'last');

for i_loop_ESU=1:length(seq_ESU)-1
    %for i_loop_ESU=1:1
    %for i_loop_ESU=101:length(seq_ESU)-1
    
    disp(['loop number ', num2str(i_loop_ESU), ' over ', num2str(length(seq_ESU)-1)])
    
    %esu indices in current loop
    loop_seq_esu=seq_ESU(i_loop_ESU):seq_ESU(i_loop_ESU+1)-1;
    %For each Esu we have corresponding ping
    i_ping_deb=EsuPingStartEnd(1,loop_seq_esu);
    i_ping_end=EsuPingStartEnd(2,loop_seq_esu);
    
    %test unicite ping dans les esus :
    diff_ping=i_ping_end-i_ping_deb;
    if sum(diff_ping) >0
        disp('one esu contain more than one ping !!!')
        disp(['Ping N ', num2str(i_ping_deb(diff_ping>0))]);
        %return
    end
    clear diff_ping
    
    %test que tous les pings sont utilisés.. sinon il faudra les supprimer
    % diff_ping_2=diff(i_ping_deb,1);
    % if sum(diff_ping_2~=1) >0
    %     disp('Some pings was skipped')
    %     ind_ping_loop_not_used= [(EsuPingStartEnd((diff_ping_2>1)+1,1)):EsuPingStartEnd(i_ping_deb(diff_ping_2>1)+1,2)-1];
    %
    % end
    % clear diff_ping_2
    
    n_esu=size(i_ping_deb,2);
    class_mask=Classif.Mask_clean(1:ind_maxdepth_EI,loop_seq_esu,:);
    % -20 dB au dessu de la première couche

    % if i_loop_ESU==1
    %     echogram_maskclean_1=ones(ind_maxsampledepth_1,max(i_ping_end)); %init to ones
    %     echogram_maskclean_2=ones(ind_maxsampledepth_2,max(i_ping_end)); %init to ones
    % else
    echogram_maskclean_1=ones(ind_maxsampledepth_1,(max(i_ping_end)-min(i_ping_deb))+1); %init to ones
    echogram_maskclean_2=ones(ind_maxsampledepth_2,(max(i_ping_end)-min(i_ping_deb))+1); %init to ones
    % end
    
    %Set to NaN shallow samples (shalower than first EI layer
    echogram_maskclean_1(1:ind_minsampledepth_1,:)=NaN;
    echogram_maskclean_2(1:ind_minsampledepth_2,:)=NaN;
    
    for i_esu=1:n_esu
        
        for i_depth=1:n_depth_EI
            
            ind_ping_loop=i_ping_deb(i_esu)-i_ping_deb(1)+1:(i_ping_end(i_esu)-i_ping_end(1)+1);
            %ind_ping_loop=i_esu;
            
            d_ind_1=(indEI_forsampl_1==i_depth);
            echogram_maskclean_1( d_ind_1 ,ind_ping_loop)=  class_mask(i_depth,i_esu,1);
            d_ind_2=(indEI_forsampl_2==i_depth);
            echogram_maskclean_2( d_ind_2 , ind_ping_loop)=  class_mask(i_depth,i_esu,2);
            %
            
            %indEI_forsampl_1(i_depth)
            %         Classif.EsuPingStartEnd(:,1)
        end
        
        
    end
    
    
    
    %gerer les pings non utilisés
    %ind_ping_loop_not_used +
    %seq_ESU(i_loop_ESU)-1
    
    %mise a jour dans echo : attention aux NaN qu'il y avait avant !!!
    Old_manual_cor_1=isnan(filt.Mask_Cleaning70(1:size(echogram_maskclean_1,1), i_ping_deb(1):i_ping_end(end)));
    echogram_maskclean_1(Old_manual_cor_1)=NaN;
    % Fish_NaN_1=isnan(echogram_maskclean_1);
    filt.Mask_Cleaning70(1:size(echogram_maskclean_1,1), i_ping_deb(1):i_ping_end(end))= echogram_maskclean_1 ;
    clear Old_manual_cor_1
    
    Old_manual_cor_2=isnan(filt.Mask_Cleaning200(1:size(echogram_maskclean_2,1), i_ping_deb(1):i_ping_end(end)));
    % Fish_NaN_2=isnan(echogram_maskclean_2);
    echogram_maskclean_2(Old_manual_cor_2)=NaN;
    filt.Mask_Cleaning200(1:size(echogram_maskclean_2,1), i_ping_deb(1):i_ping_end(end))= echogram_maskclean_2 ;
    clear Old_manual_cor_2
    
    %mettre le filtering à NaN pour les echantillons sous la profondeur max de
    %l'Ei
    % filt.Mask_Cleaning70(ind_maxsampledepth_1+1:end, i_ping_deb(1):i_ping_end(end))=NaN;
    % filt.Mask_Cleaning200(ind_maxsampledepth_2+1:end, i_ping_deb(1):i_ping_end(end))=NaN;
    
    
    for i_sampledepth=ind_maxsampledepth_1+1:depth_samples_1
        filt.Mask_Cleaning70(i_sampledepth, i_ping_deb(1):i_ping_end(end))=NaN;
        
    end
    for i_sampledepth=ind_maxsampledepth_2+1:depth_samples_2
        filt.Mask_Cleaning200(i_sampledepth, i_ping_deb(1):i_ping_end(end))=NaN;
    end
    %mettre le filtering à 1 pour les echantillons juste au dessus de la
    %premiere couche. Comme ca la production de clean hac dans matecho ne
    %supprimera pas les ping n'ayant aucun echantillon
    filt.Mask_Cleaning70(ind_sample_20dB_1, i_ping_deb(1):i_ping_end(end))=1;
    filt.Mask_Cleaning200(ind_sample_20dB_2, i_ping_deb(1):i_ping_end(end))=1;
    
    disp(['OK from ping ', num2str(i_ping_deb(1)), ' to ', num2str(i_ping_end(end)),]);
    
    %input('tape une touche')%
    
end % for i_loop_ESU=1:length(seq_ESU)

%mettre -20dB dans l'echogram juste au dessus de la
%premiere couche. Comme ca la production de clean hac dans matecho ne
%supprimera pas les ping n'ayant aucun echantillon
% echo.Echogram70(ind_sample_20dB_1,:)=-20;
% echo.Echogram200(ind_sample_20dB_2,:)=-20;


