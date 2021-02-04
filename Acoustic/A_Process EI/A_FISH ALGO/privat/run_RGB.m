%build dB limits
lim01=[Low_Limits(1), Low_Limits(1)+Ranges(1)];
lim02=[Low_Limits(2), Low_Limits(2)+Ranges(2)];
lim03=[Low_Limits(3), Low_Limits(1)+Ranges(3)];

%Initialise popup menu to first frequency (Red)
Frequ_work=1;
depth=depth_surface(1,1:id_maxD);
Ping_index=[1:1:length(sv01(:,1))];
%call RGB_gui function
    RGB_gui(Ping_index, depth, lim01, lim02, lim03, sv01, sv02, sv03, Frequ_work,Freq,Night_Day_sel,fig_number)

