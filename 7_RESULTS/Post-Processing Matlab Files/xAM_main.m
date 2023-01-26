%% A.M. 23-02-2020 - This script reads and displays the results1 contained in the HDF5 data
% clear all
% clc
for slicenum1 = [2]%[12,3,4,5,6]
   
        clearvars -except slicenum1

        %% Initialize parameters for all the significant variables
        [medium,domain,dslice,file,plin,pnl] = xAM_Init(slicenum1);
        %%
    %% Bubble Parameters
    Bubble.mindist           = 10e-6;
    Lin=Bubble;
    [Bubble] = xAM_BubbleCluster_LocCon(medium,domain,file,Bubble,'Bubble');
    [Lin] = xAM_BubbleCluster_LocCon(medium,domain,file,Lin,'Lin');
    %%  Bubble Cluster Positions
    [domain] = xAM_BubbleCluster_LocPlot(domain,dslice,file,Bubble,'Microbubble','.',8);
    %%
    [domain] = xAM_BubbleCluster_LocPlot(domain,dslice,file,Lin,'Microbubble','.',10);
    %% PlayMovies
    xAM_PlayMovies(medium,domain,dslice,file,plin,pnl,Bubble,Lin)
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
    disp('Bubble Parameters ...' ) ; disp(Bubble);
end