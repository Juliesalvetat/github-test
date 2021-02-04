%build dB limits
lim01=[Low_Limits(1), Low_Limits(1)+Ranges(1)];
lim02=[Low_Limits(2), Low_Limits(2)+Ranges(2)];
lim03=[Low_Limits(3), Low_Limits(1)+Ranges(3)];

%Initialise popup menu to first frequency (Red)
Frequ_work=1;
depth=depth_surface(1,1:id_maxD);
Ping_index=[1:1:length(sv01(:,1))];
%call RGB_gui function
%    RGB_gui(Ping_index, depth, lim01, lim02, lim03, sv01, sv02, sv03, Frequ_work,Freq,Night1Sunrise2Day3Sunset4,fig_number)


%%%%%%%%  Transform dB in color indices %%%%%%%%%%
fit01=polyfit(lim01,[0,255],1);
color01=sv01*fit01(1)+fit01(2);

fit02=polyfit(lim02,[0,255],1);
color02=sv02*fit02(1)+fit02(2);

fit03=polyfit(lim03,[0,255],1);
color03=sv03*fit03(1)+fit03(2);

%%%%% Build a image with indices from 0 to 255  %%%%%%%%%%
echocolor=cat(3,uint8(round(color01')),uint8(round(color02')),uint8(round(color03')));
echocolor(isnan(echocolor))=0;
echocolor(echocolor<0)=0;
echocolor(echocolor>255)=255;

 
 %%%%%% Figure 1 is echogram in RGB color scale %%%%%%%%%%%%%%
f1=figure(1);
set(f1,'Position',[0 20 1800 700],'Name','Three frequencies echogram in RGB color scale"');




% image(Ping_index, depth , echocolor,...
% 'CDataMapping','scaled');
positionVector1 = [0.1, 0.15, 0.72, 0.8];
sp1=subplot('Position',positionVector1);

imagesc(Ping_index, depth , echocolor);
set(gca,'xtick',[]); 
ColorNightDay=[0 0 0;1 1 0;1 1 1;1 0 0]; 
title(['R : ', num2str(Freq(1)), ' kHz, G : ', num2str(Freq(2)), ' kHz, B : ', num2str(Freq(3)), ' kHz'])
xlabel('ESU')
ylabel('Depth (m)');
