buton=1;
esu=[];
ind_depth=[];

continue_selec=1;
id_polygon=1;

while continue_selec
    
    clear x y
    
    
    while buton==1
        
        [x,y,buton]=ginput(1);
        esu=[esu,x];
        ind_depth=[ind_depth,y];
        hold on
        if length(esu)>=2
            plot([esu(end-1),esu(end)],[ind_depth(end-1),ind_depth(end)],'r');
        end
        
        
    end
    hold off
    
    polygon(id_polygon).esu=esu;
    polygon(id_polygon).ind_depth=ind_depth;
    esu=[];
    ind_depth=[];

    %ASk user to continue or not
    choice = questdlg('Would you like to do an other selection ?', ...
        'Yes','Yes','No thank you','No thank you');
    % Handle response
    switch choice
        case 'Yes'
            continue_selec=1;
            id_polygon=id_polygon+1;
            buton=1; hold off ;
            
            %set(gcf,'HitTest','off');
            str='No';
            while ~strcmp(str,'')
            str = input('Type Enter when ready after zooming/unzooming on figure)','s');
            end
            
            
        case 'No thank you'
            continue_selec=0;
    end %end while buton==1
    
    
end %end while continue_selec