function RGB_gui(Ping_index, depth, lim01, lim02, lim03, sv01, sv02, sv03, Frequ_work,Freq,Night1Sunrise2Day3Sunset4,fignum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read user parameters for processing
%D_echocolors_Define_param
%we need to read it again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
f1=figure(fignum);
set(f1,'Position',[0 20 1800 700],'Name','Three frequencies echogram in RGB color scale"');

positionVector00 = [0.9, 0.1, 0.05, 0.5];
sp00=subplot('Position',positionVector00);

size_mark=150;
hold on
for i_dyn=lim01(1):1:lim01(end) %pour chaque dB de la dynamique de couleur 1  
  scatter(1,i_dyn,size_mark,[(round(i_dyn*fit01(1)+fit01(2)))/255,0, 0],'s','fill');
end
for i_dyn=lim02(1):1:lim02(end) %pour chaque dB de la dynamique de couleur 2  
  scatter(2,i_dyn,size_mark,[0,(round(i_dyn*fit01(1)+fit01(2)))/255, 0],'s','fill');
end
for i_dyn=lim01(1):1:lim01(end) %pour chaque dB de la dynamique de couleur 3
  scatter(3,i_dyn,size_mark,[0,0,(round(i_dyn*fit01(1)+fit01(2)))/255],'s','fill');
end
for i_dyn=max([lim01(1),lim02(1),lim03(1)]):1:min([lim01(end),lim02(end),lim03(end)]) % Les trois couleurs
    scatter(4,i_dyn,size_mark,[(round(i_dyn*fit01(1)+fit01(2)))/255, (round(i_dyn*fit01(1)+fit01(2)))/255 , (round(i_dyn*fit01(1)+fit01(2)))/255],'s','fill');
end

xlim([0.5 4.5])
ylim([min([lim01(1),lim02(1),lim03(1)])-5;  max([lim01(end),lim02(end),lim03(end)])+5])
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'R', 'G', 'B', 'RGB'});
xlabel('Color')
ylabel('dB')
hold off
set(f1,'Visible','on');
set(sp00,'NextPlot','replace');


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
positionVector2 = [0.1, 0.04, 0.72, 0.05];
sp0=subplot('Position',positionVector2);


imagesc(Night1Sunrise2Day3Sunset4); caxis([1 4]); 
colormap(sp0,ColorNightDay); 
set(gca,'xtick',[]); set(gca,'ytick',[]); 
xlabel('Night->Black / Sunrise->Yellow / Day->White / Sunset->Red'); 
             
           

%%%%%%%%%%% Build some uicontrol %%%%%%%%%%%%%%%%%%%%%%ù

hp1 = uipanel('Title','Color adjust','FontSize',11,...
             'BackgroundColor','white',...
             'Position',[.87 .67 .12 .25]);

htext = uicontrol('Style','text','String','Color (Frequency)','Units','normalized',...
           'Position',[.88,.79,.1,.1]);
       
 hpopup = uicontrol('Style','popupmenu',...
           'String',{'Red','Green','Blue','All'},'Value',Frequ_work,'Units','normalized',...
           'Position',[.88,.75,.1,.1]);
      
hRplus    = uicontrol('Style','pushbutton',...
             'String','Color +','Units','normalized',...
             'Position',[.88,.73,.1,.03],'Callback',{@Cplus});
hRmoins    = uicontrol('Style','pushbutton',...
             'String','Color -','Units','normalized',...
             'Position',[.88,.68,.1,.03],'Callback',{@Cmoins});

         




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a color histogram %%%%%%%%%
% f3=figure(3);
% set(f3,'Position',[875 150 390 750],'Name','Color histogram"');
% B=echocolor(:,:,3);
% V=echocolor(:,:,2);
% R=echocolor(:,:,1);
% hist([B(:),V(:), R(:)],[0:2:255])
% [n,xout]=hist([B(:),V(:), R(:)],[0:2:255]);
% maxY=max((ceil(max(n(5:123,:))/8))*10);
% ylim([0 maxY])
% xlim([0 255])
% title(['Red : ', num2str(Freq(1)), ' kHz, Green : ', num2str(Freq(2)), ' kHz, Blue : ', num2str(Freq(3)), ' kHz'])
% xlabel('Color index [0-255]')
% ylabel('Histogram');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
 %%%%%% Function to increase color %%%%%%%%%%%%%  
 function Cplus(source,eventdata)
     Frequ_work=get(hpopup ,'Value');    %read frequency (color) from popupmenu      
     switch Frequ_work
         case 1 % 2dB less for limits of first frequency
            RGB_gui(Ping_index, depth, lim01-3, lim02, lim03, sv01, sv02, sv03,1,Freq,Night1Sunrise2Day3Sunset4,fignum)
         case 2 % 2dB less for limits of 2nd frequency
            RGB_gui(Ping_index, depth, lim01, lim02-3, lim03, sv01, sv02, sv03,2,Freq,Night1Sunrise2Day3Sunset4,fignum)
         case 3 % 2dB less for limits of 3rd frequency
            RGB_gui(Ping_index, depth, lim01, lim02, lim03-3, sv01, sv02, sv03,3,Freq,Night1Sunrise2Day3Sunset4,fignum)
         case 4 %2dB less for limits at all frequencies
            RGB_gui(Ping_index, depth, lim01-3, lim02-3, lim03-3, sv01, sv02, sv03,4,Freq,Night1Sunrise2Day3Sunset4,fignum)
     end
     disp(' ')
     disp('Updating data. Please wait.... ')
     disp(' ')
    
 end

 %%%%%% Function to decrease color %%%%%%%%%%%%%  
 function Cmoins(source,eventdata)
     Frequ_work=get(hpopup ,'Value'); %read frequency (color) from popupmenu      
     switch Frequ_work
         case 1 %add 2dB in limits of first frequency
            RGB_gui(Ping_index, depth, lim01+3, lim02, lim03, sv01, sv02, sv03,1,Freq,Night1Sunrise2Day3Sunset4,fignum)
         case 2 %add 2dB in limits of 2nd frequency
            RGB_gui(Ping_index, depth, lim01, lim02+3, lim03, sv01, sv02, sv03,2,Freq,Night1Sunrise2Day3Sunset4,fignum)
         case 3 %add 2dB in limits of 3rd frequency
            RGB_gui(Ping_index, depth, lim01, lim02, lim03+3, sv01, sv02, sv03,3,Freq,Night1Sunrise2Day3Sunset4,fignum)
         case 4 %add 2dB in limits at all frequencies
            RGB_gui(Ping_index, depth, lim01+3, lim02+3, lim03+3, sv01, sv02, sv03,4,Freq,Night1Sunrise2Day3Sunset4,fignum)
     end
     disp(' ')
     disp('Updating data. Please wait.... ')
     disp(' ')
    
 end


end %end of function RVB_gui
       
         


