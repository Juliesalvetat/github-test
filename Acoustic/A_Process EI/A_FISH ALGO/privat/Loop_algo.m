close all;
ek5=ek500();

Sv_surface=echo.Sv_surface(1:id_maxD,ESU_SELection,:);
Mask_clean=echo.Mask_clean(1:id_maxD,ESU_SELection,:);
%Cleaning mask from manual EI corrections
Sv_surface(isnan(Mask_clean))=NaN;
%Seuillage des Sv. L'algo ne se fait qu'avec au moins un des 2 Sv superieur à ce seuil
mask_LowSv_algo=(Sv_surface(:,:,1)<Sv_LowThresh) & (Sv_surface(:,:,2)<Sv_LowThresh);
mask_LowSv_algo=repmat(mask_LowSv_algo,[1 1 length(FrequencySort)]);
Sv_surface(mask_LowSv_algo)=NaN;

sv01=Sv_surface(:,:,idF01)';
sv02=Sv_surface(:,:,idF02)';
sv03=Sv_surface(:,:,idF03)'; %ONLY FOR RGB PLOts



Night_Day_sel=Night1Sunrise2Day3Sunset4(1,ESU_SELection);
depth_surface=depth_surface(1,1:id_maxD);

%plot option  for sv01
% figure()
% pcolor(sv01');shading flat; set(gca,'Clim',[-80 -40],'YDir','reverse');colorbar;

%PLOT RGB before calssification
if Process_stop_loop
    Bottom=echo.Bottom;
    %     fig_number=98;
    %     run_RGB
    
    ind_depth=find(depth_surface<max(Bottom(ESU_SELection)+5));
    
    figure(1)
    set(gcf,'position',get(0,'screensize'))
    subplot(2,1,1)
    imagesc(size(sv01,1),depth_surface(ind_depth),sv01(:,ind_depth)');shading flat; set(gca,'YDir','reverse','CLim',[-80 -40]);colorbar;colormap(ek5);
    title('Original 70 kHz')
    subplot(2,1,2)
    imagesc(size(sv01,1),depth_surface(ind_depth),sv02(:,ind_depth)');shading flat; set(gca,'YDir','reverse','CLim',[-80 -40]);colorbar;
    title('Original 200 kHz')
end

masknan01=isnan(sv01);
masknan02=isnan(sv02);

%nombre ech pour calculer la variance
Nsamp_dim3=(Size_filt_vertical*2+1)*(Size_filt_horiz*2+1);

%Build a 3D matrix to compute variance of each Sv :
%cette  variable est "élargie" de la taille du filtre
Decal_Sv01=NaN([size(sv01,1)+(2*Size_filt_horiz),size(sv01,2)+(2*Size_filt_vertical),Nsamp_dim3]);
Decal_Sv02=NaN([size(sv02,1)+(2*Size_filt_horiz),size(sv02,2)+(2*Size_filt_vertical),Nsamp_dim3]);

Ech_central_vert=Size_filt_vertical+1;
Ech_central_hor=Size_filt_horiz+1;
isamp_3D=0;

%pour chaque echantillon du filtre on va concatener dans la dimension 3
for isampH=-Size_filt_horiz:Size_filt_horiz
    Hdecal=Ech_central_hor+isampH;
    
    for isampV=-Size_filt_vertical:Size_filt_vertical
        isamp_3D=isamp_3D+1;
        Vdecal=Ech_central_vert+isampV;
        Decal_Sv01(Hdecal:size(sv01,1)+Hdecal-1,Vdecal:size(sv01,2)+Vdecal-1,isamp_3D)=sv01;
        Decal_Sv02(Hdecal:size(sv02,1)+Hdecal-1,Vdecal:size(sv02,2)+Vdecal-1,isamp_3D)=sv02;
    end
    
end %isamp=1:Nsamp_dim3



Var_Sv01=var(Decal_Sv01,1,3,'omitnan');
Var_Sv02=var(Decal_Sv02,1,3,'omitnan');
%if Matlab version is below V2016, use next lines
% Var_Sv01=nanvar(Decal_Sv01,3);
% Var_Sv02=nanvar(Decal_Sv02,3);

clear Decal_Sv01 Decal_Sv02


%delete rows ans columns that have been created before
Var_Sv01(1:Size_filt_horiz,:)=[];
Var_Sv01((end-Size_filt_horiz+1):end,:)=[];
Var_Sv02(1:Size_filt_horiz,:)=[];
Var_Sv02((end-Size_filt_horiz+1):end,:)=[];

Var_Sv01(:,1:Size_filt_vertical)=[];
Var_Sv01(:,(end-Size_filt_vertical+1):end)=[];
Var_Sv02(:,1:Size_filt_vertical)=[];
Var_Sv02(:,(end-Size_filt_vertical+1):end)=[];


%Apply original mask of NaN
Var_Sv01(masknan01)=NaN;
Var_Sv02(masknan02)=NaN;

% If necessary, plot varianc of each Sv
if Process_stop_loop
    figure(2)
    set(gcf,'position',get(0,'screensize'))
    subplot(2,1,1)
    imagesc(size(sv01,1),depth_surface(ind_depth),Var_Sv01(:,ind_depth)');shading flat; set(gca,'YDir','reverse','CLim',[0 150]);colorbar;colormap(jet);
    %imagesc(Var_Sv01');shading flat; set(gca,'YDir','reverse','CLim',[0 4]);colorbar;
    subplot(2,1,2)
    imagesc(size(sv01,1),depth_surface(ind_depth),Var_Sv02(:,ind_depth)');shading flat; set(gca,'YDir','reverse','CLim',[0 150]);colorbar;colormap(jet);
    %imagesc(Var_Sv02');shading flat; set(gca,'YDir','reverse','CLim',[0 4]);colorbar;
end



%avec Sv01 et Sv02 : seuil sur la la variance
mask_fish_var=Var_Sv01>thresh_var_1 | Var_Sv02>thresh_var_2;
%mask_fish_var= Var_Sv02>thresh_var_2;

%Option to plot mask on variance
% figure()
% pcolor(double(mask_fish_var'));shading flat; set(gca,'YDir','reverse');colorbar;



%enlever les echos forts de nofish et les rajouter au masque des fortes
%variances
mask_fish= ((mask_fish_var & sv01 > Tresh_sv1) | (mask_fish_var & sv02 > Tresh_sv2))| sv01+sv02 > Thresh_sum;

%mask no vish = inverse mask fish
mask_no_fish=~mask_fish;

% Save matrice before applying NoFish mask
Save_Sv01=sv01;
Save_Sv02=sv02;
Save_Sv03=sv03;

%NO FISH DATA%%%%%%%%%%%%%%%%%%%%%
sv01(~mask_no_fish)=NaN;
sv02(~mask_no_fish)=NaN;
sv03=sv02;
%RGB for NO Fish
if Process_stop_loop
    %fig_number=99;
    %run_RGB
    %     figure(3)
    %     subplot(2,1,1)
    %     Sv_noFish=Sv_surface;
    %     Sv_noFish(mask_fish')=NaN;
    %     imagesc(Sv_noFish(ind_depth,ESU_SELection,1));shading flat; set(gca,'YDir','reverse','CLim',[-80 -40]);colorbar;colormap(ek5);
    %     %ylim([0 depth_surface(ind_depth(end))])
    %     title('No Fish algorithm result - 70 kHz ')
end

%%%%%% FISH DATA
%
sv01=Save_Sv01; sv02=Save_Sv02; sv03=Save_Sv03;
clear Save_Sv01 Save_Sv02 Save_Sv03
sv01(~mask_fish)=NaN;
sv02(~mask_fish)=NaN;
sv03=sv02;
%RGB for Fish
if Process_stop_loop
    %fig_number=100;
    %run_RGB
    figure(3)
    set(gcf,'position',get(0,'screensize'))
    subplot(2,1,1)
    imagesc(size(sv01,1),depth_surface(ind_depth),sv01(:,ind_depth)');shading flat; set(gca,'YDir','reverse','CLim',[-80 -40]);colorbar;colormap(ek5);
    title('Fish algorithm result - 70 kHz ')
    subplot(2,1,2)
    imagesc(size(sv01,1),depth_surface(ind_depth),sv02(:,ind_depth)');shading flat; set(gca,'YDir','reverse','CLim',[-80 -40]);colorbar;colormap(ek5);
    title('Fish algorithm result - 200 kHz ')
    
    % ylim([0 ind_depth(end)])
    
end

