%function to save and print data for a group
    function save_print(mask,nfigure)
        defaultfilename=['Group', num2str(nfigure-3),'.mat'];
        [FileNameGroup,PathNameGroup] = uiputfile('*.mat','Enter file name to save group mask',defaultfilename) ;
        save([PathNameGroup,FileNameGroup],'mask', '-mat');
        
        disp(' ')
        disp('Mask saved in : ')
        disp([PathNameGroup,FileNameGroup])
        disp(' ')
        
        figure(nfigure)
        subplot(2,1,1)
        xlabel('') %suppress label before printing ('ESU' overwrited with text)
        FileNameFig=[FileNameGroup(1:end-3), 'pdf'];
        print([PathNameGroup,FileNameFig],'-dpdf');
        xlabel('ESU') 


        disp(' ')
        disp('Figure saved in : ')
        disp([PathNameGroup,FileNameFig])
        disp(' ')

    end
