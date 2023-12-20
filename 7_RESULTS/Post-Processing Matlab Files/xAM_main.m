%% A.M. 23-02-2020 - This script reads and displays the results1 contained in the HDF5 data
% clear all
% clc
for slicenum1 = [2]%[12,3,4,5,6]
   
        clearvars -except slicenum1 q

        %% Initialize parameters for all the significant variables
        [medium,domain,dslice,file,plin,pnl] = xAM_Init(slicenum1);
        %%
    %% Bubble Parameters
    Bubble.mindist           = 10e-6;
    for ipop = 1:4
        name_pop={'Microbubble','Linear Scatterer','Microbubble','Linear Scatterer'};
        [ScPop{ipop}] = xAM_BubbleCluster_LocCon(medium,domain,file,Bubble,'Bubble',ipop);
        %  Bubble Cluster Positions
        [domain] = xAM_BubbleCluster_LocPlot(domain,dslice,file,ScPop{ipop},name_pop{ipop},'.',8);
        clear name_pop
    end
    %% PlayMovies
    xAM_PlayMovies(medium,domain,dslice,file,plin,pnl,ScPop)
    %% CREATING MATRICES FOR NONLINEAR AND CONTRAST DATA
    [pnl,plin] = xAM_ConPlot(medium,domain,dslice,file,plin,pnl);
    %%
    fclose('all');  
    disp(' ')
    disp('Medium Parameters ...' ) ; disp(medium)  ;
    disp('Domain Parameters ...' ) ; disp(domain ) ;
    disp('Slice Parameters ...' ) ; disp(dslice ) ;
    disp('File Parameters ...' ) ; disp(file )  ;
    disp('Linear Pressure Parameters ...' ) ; disp(plin ) ;
    disp('Nonlinear Pressure Parameters ...' ) ; disp(pnl )   ;
    disp('Bubble Parameters ...' ) ; disp(ScPop);
end