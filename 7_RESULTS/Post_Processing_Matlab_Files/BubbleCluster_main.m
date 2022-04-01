%% A.M. 23-02-2020 - This script reads and displays the results1 contained in the HDF5 data
% clear all
% fclose('all')
% close all;
for slicenum1 = [2]%[12,3,4,5,6]
    
    clearvars -except slicenum1 atten atten_fund Loc pnl7
    %     close all
    clc
    %% Initialize parameters for all the significant variables
    [medium,domain,dslice,file,plin,pnl,Bubble] = BubbleCluster_Init(slicenum1);
    %%
%     %% Bubble Location and Pressure
    [Bubble] = BubbleCluster_LocCon(medium,domain,file,Bubble,'Bubble');
    PS=Bubble ; domainPS = domain; domainPS.beamiterations=0;[PS] = BubbleCluster_LocCon(medium,domainPS,file,PS,'PS');
    %%  Bubble Cluster Positions
    [domain] = BubbleCluster_LocPlot(domain,dslice,file,Bubble,'Microbubble','o',3);
    [domain] = BubbleCluster_LocPlot(domain,dslice,file,PS, 'Point Source','kx',0.8);
    
    %% PlayMovies
     Bubble_PlayMovies(medium,domain,dslice,file,plin,pnl,Bubble)
    %% CREATING MATRICES FOR NONLINEAR AND CONTRAST DATA
    [pnl,plin] = BubbleCluster_ConPlot(medium,domain,dslice,file,plin,pnl,Bubble);
%     atten(dslice.num,:) = pnl.attenuation(dslice.num,:);
%     atten_fund(dslice.num,:) = pnl.attenuation_fund(dslice.num,:);
    %% Compare Analytic solution with simulation results
    BubbleCluster_Compare(medium,domain,dslice,file,plin,pnl,Bubble);
    %% Plot residuals
    Bubble_ResidualsPlot(medium,domain,dslice,file,plin,pnl,Bubble,[1E2 1E3 1E4],'../test/Phased_2E5Pa_2micron_1MHz_')
    
    fclose('all');
    %%
    
    disp(' ')
    disp('Medium Parameters ...' ) ; disp(medium)  ;
    disp('Domain Parameters ...' ) ; disp(domain ) ;
    disp('Slice Parameters ...' ) ; disp(dslice ) ;
    disp('File Parameters ...' ) ; disp(file )  ;
    disp('Linear Pressure Parameters ...' ) ; disp(plin ) ;
    disp('Nonlinear Pressure Parameters ...' ) ; disp(pnl )   ;
    disp('Bubble Parameters ...' ) ; disp(Bubble);
    
    %%
    
    if strcmp(file.plot_attenslices, 'yes')
        f9 = figure('WindowState','maximized');
        slice= 6;
        hold on;
        plot(atten(slice:end,1),abs(atten(slice:end,2)),'x--k',atten_fund(slice:end,1),abs(atten_fund(slice:end,2)),'o--k')
        %         plot(atten9(5:end-1,1),abs(atten9(5:end-1,2)),'x:k',atten_fund9(5:end-1,1),abs(atten_fund9(5:end-1,2)),'o:k')
        
        xline(min(Bubble.LocRange(3,:)),'Color','blue');
        xline(max(Bubble.LocRange(3,:)),'Color','blue');
        xlim([Bubble.LocRange(3,1) Bubble.LocRange(3,2)])
        xlabel('Bubble Cluster Position [mm]')
        ylabel('Attenuation [dB/cm]')
        title(' Attenuation of the wave propagated through bubble cluster')
        legend('Full Spectrum','Fundamental','Location','southeast')
        %         legend('Full Spectrum, 6MHz','Fundamental, 6MHz','Full Spectrum, 9MHz','Fundamental, 9MHz')
        grid on;
        set(gca,'FontSize',25)
        xlim([min(Bubble.LocRange(3,:))-1 max(atten(slice:end,1))+1])
    end
end