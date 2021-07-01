% ***********************************************************************************************************
% Compare Analytic solution with simulation results
%
% ***************************************************************
% BubbleCluster_Compare(medium,domain,slice,file,plin,pnl,Bubble)
% ***************************************************************
% Compare the results generated from the matlab ODE solver with the ODEPACK Fortran90 solver
% for the Marmottant model , which describes how a micobubble behaves when excited by an acoustic wave
% Compare the results of the simulation with the analytic expression of the scattered pressure of one
% scatterer.
%
% Inputs:
% medium        struct      Contains all the medium properties, e.g Density, Speed of Sound etc
% domain        struct      Contains all the computational domain properties, e.g space and time step
% slice         struct      Contains all the information about the slice ,e.g Dimension, Position, Index etc
% file          struct      Contains all the loaded file properties ,e.g directory and root name etc
% plin          struct      Contains all the pressure values of the plin slice
% pnl           struct      Contains all the pressure values of the pnl slice
% Bubbles       struct      Contains all properties of the Scatterers' Cluster
%
% ******************
% AM, version 200723
% ******************
%
% ***********************************************************************************************************


function BubbleCluster_Compare(medium,domain,dslice,file,plin,pnl,Bubble)
disp(' ')
disp('Calculating the difference between analytic and simulation solution...')
if (Bubble.N <=1E5 && Bubble.N>0)
    
    Sim_OSFactor=1;
    %     % Load pressure values from file
    %             filename_pres = split(ls([file.dirname '/Bubbles/' 'output_file_pressure_*.txt']));
    %             filename_pres = filename_pres(1:end-1);
    %             for ifile = 1:size(filename_pres)
    %                 fileID_pres = fopen(filename_pres{ifile},'r');
    %                 size_pres = [1 Inf];
    %                 formatSpec = '%f';
    %                 BubbleVal_loaded = fscanf(fileID_pres,formatSpec,size_pres);
    %
    %                 n_pad = 4;
    %         %         BubbleVal = BubbleVal_loaded(1+size(BubbleVal_loaded,1)/2*(iBubble-1):end*iBubble/2);
    %                 BubbleVal = BubbleVal_loaded(2:(1+n_pad)*2*domain.tdimpar+1);
    %
    %                 BubbleTime = BubbleVal(1:domain.tdimpar)-domain.dtpar;
    %                 BubblePres = BubbleVal(domain.tdimpar + 1 : 2*domain.tdimpar);
    %                 BubbleTimePad = BubbleVal(2*domain.tdimpar + 1:(n_pad+2)*domain.tdimpar)-domain.dtpar;
    %                 BubblePresPad = BubbleVal((n_pad+2)*domain.tdimpar + 1 : end);
    %                 figure; plot(BubbleTime,BubblePres);hold on; plot(BubbleTimePad,BubblePresPad)
    %
    %                 % Create Simulation Contrast, correct if n_padded in INCS code
    %                 n_pad = length(BubblePres)/domain.tdimpar;
    %                 BubbleTimeSim = interp(domain.tpar,n_pad);
    %                 Sim_OSFactor = n_pad*length(Bubble.Contrast)/length(BubbleTimeSim);
    %                 if Sim_OSFactor>1 ; BubbleCon_Sim = downsample(Bubble.Contrast,Sim_OSFactor/n_pad); else ; BubbleCon_Sim = interp(Bubble.Contrast,n_pad/Sim_OSFactor); end
    %             end
    %% ============== Comparison of time signature based on ODE Solvers ====================================
    if (strcmp(file.scatterer,'active'))
        
        V_dd_norm = Bubble_SimAnalytic(BubblePresPad,BubbleTimePad,medium.c0,medium.rho0,medium.freq0);         % Solution of RP Solver with ode45 MATLAB solver
        BubbleCon_an = V_dd_norm*medium.rho0;                                                             % Contrast Source Term
        %========================== Save plot ===================================
        if (strcmp(file.saveplot,'yes'))
            saveas(gcf,[file.savedir,'/',file.dirname,'_bubble_motion_slice',int2str(dslice.num),'.jpg'])
        end
        
        f5 = figure('WindowState','maximized');
        ax = gca; hold on;
        plot(BubbleTimeSim,BubbleCon_Sim,'--','color','k');
        plot(BubbleTimePad,BubbleCon_an,'-','color',[17 17 17]./255);
        legend('Simulation', 'Analytic')
        ax.YLabel.String = 'Pressure [Pa]';
        BubbleCon_Sim(isnan(BubbleCon_Sim)) = 0 ;
        %         RRMS_ErrorG = sqrt( sum(sum((BubbleCon_Sim - BubbleCon_an).^2))./sum(BubbleCon_an.^2));
        set(ax,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [microsec]';
    end
    %     end
    %% =============================================== Compare simulation with analutical Results =======================================================
    % This is working for 1 or 2 Bubbles and up to 1st order multiple scattering
    % it can easily expanded to more but these are enough for the validation of the code
    %     close all
    if (Bubble.N<=1E5 )
        if (Bubble.N==2); domain.beamiterations=1;end
        if (domain.beamiterations>2) ;domain.beamiterations=2;end
        
        OSFactor = 	1;
        
        n_pad_fft = [];%2^nextpow2(domain.tdimpar);
        
        k_init      =   1;
        k_end       =   domain.dimlen(1);
        m_init      =   1;%domain.dimlen(3);
        m_end       =   domain.dimlen(3);
        %         % %
        %                 k_init      =   58;
        %                 k_end       =   k_init;
        %                 m_init      =   77;
        %                 m_end       =   m_init;
        RRMS_ErrorG_semi_an = zeros(k_end,m_end);
        RRMS_ErrorG_an = zeros(k_end,m_end);
        
        %======== Driving Frequency ======================
        freq        = medium.freq0;
        w_driving   = 2*pi*freq;
        tpar_N      = interp(domain.tpar,OSFactor);
        Inc_Press   = 1E4*exp( - ((tpar_N-6/freq)/(3/freq/2)).^2).* sin(w_driving*(tpar_N-6/freq)).*(1+sign(tpar_N))/2;
        Inc_Press   = repmat(Inc_Press,Bubble.N,1);
        omega       = 2*pi/(domain.tdimpar*domain.dtpar)*((0:domain.tdimpar*OSFactor-1)-ceil((domain.tdimpar*OSFactor-1)/2));
        omega       = [omega zeros(1,n_pad_fft-domain.tdimpar)];
        omega       = repmat(omega,Bubble.N,1);
        
        % Make denser the matrix to increase accuracy
        pdiff  = pnl.contrastdata;  % If you do not save contrastdata,  pdiff = plin.data - pnl.data
        
        %======= Initialize Arrays ========================
        pnl.simN                = zeros(domain.tdimpar*OSFactor,1);
        pnl.anN                 = zeros(domain.tdimpar*OSFactor,Bubble.N);
        pnl.semi_anN            = zeros(domain.tdimpar*OSFactor,Bubble.N);
        
        % Iterate over k in x axis and m in z axis. This implementation is done for a y slice so the index in y equals always to the closest point of the slice position (SliceIndex(2))
        for k = k_init: k_end
            for m = m_init:m_end
                P_scattered = zeros(Bubble.N,domain.tdimpar*OSFactor);
                % Iterate to include higher orders of multiple scattering
                for i = domain.beamiterations: domain.beamiterations
                    %======= Initialize Arrays ========================
                    pnl.simN_cluster        = zeros(domain.tdimpar*OSFactor,1);
                    pnl.anN_cluster         = zeros(domain.tdimpar*OSFactor,1);
                    pnl.semi_anN_cluster    = zeros(domain.tdimpar*OSFactor,1);
                    
                    % If 0 order then no scattering between scatterers, else include scattering in the for loop
                    % this is done only for 2 bubbles, otherwise another more complicated idea should be implemented
                    if (i==1); Bubbles = [1:Bubble.N]; else ;Bubbles = [Bubble.N:-1:1];end
                    for j = 1: i
                        % Initialize for higher order multiple scattering . This is done because on previous j's scattered pressure is calculated
                        if (j>1) ;Bubbles = [1:Bubble.N] ;pnl.anN_cluster = zeros(domain.tdimpar*OSFactor,1); pnl.semi_anN_cluster=pnl.anN_cluster;end
                        % Find the position of interest which is the position of slice
                        % Find the grid points with which the bubble has the smallest distance
                        %                         [~,SliceIndex(1)] = min(abs(Bubble.LocGlob(Bubbles(iBubble),domain.dimval(2)) - domain.par{domain.dimval(2)}));
                        [~,SliceIndex(2)] = min(abs(dslice.pos(dslice.num)-domain.par{domain.dimval(1)}));
                        if (i==1 || j>1) ;SliceIndex(1) = k ;  SliceIndex(3) = m; else; pnl.anN_cluster = zeros(domain.tdimpar*OSFactor,1); pnl.semi_anN_cluster=pnl.anN;  end
                        
                        %================================================== SIMULATION RESULTS =============================================================
                        
                        % The check for the time signature for microbubbles is done from the previous step of the comparison of the ode solvers
                        pnl.sim = squeeze(pdiff(:,SliceIndex(1),SliceIndex(3)));
                        pnl.simN_cluster = pnl.sim';%interp(pnl.sim,OSFactor)';
                        %================================================== GLOBAL VARIABLES =============================================================
                        % This is used to show the propagation in the horizontal axis,This is because of the results of INCS
                        r_from_scatterer = sqrt((Bubble.LocGlob(:,domain.dimval(2)) - domain.par{domain.dimval(2)}(SliceIndex(1))).^2+(Bubble.LocGlob(:,domain.dimval(1)) - dslice.pos(dslice.num)).^2 ...
                            +(Bubble.LocGlob(:,domain.dimval(3))-domain.par{domain.dimval(3)}(SliceIndex(3))).^2)*1e-3;
                        z_shift = (Bubble.LocGlob(:,domain.dimval(3))-domain.par{domain.dimval(3)}(SliceIndex(3)) )*1e-3;
                        
                        %                         if (i>1 && j==1); r_from_scatterer
                        %                         z_shift = (Bubble.LocGlob(iBubble,domain.dimval(3))-Bubble.LocGlob(Bubbles(iBubble),domain.dimval(3)) )*1e-3; end
                        
                        phase_shift = exp(-1i*omega.*z_shift/medium.c0);
                        Greens_function = (exp(-1i*omega.*r_from_scatterer/medium.c0)./(4*pi*r_from_scatterer)); % 1/(4πr) and the time shift due to wave travel
                        %================================================== SEMI-ANALYTIC RESULTS SHIFTED =============================================================
                        %                         pnl.semi_an = Bubble.Contrast ; % Semi - analytic because checking only the Green's function
                        %                         pnl.semi_an = real(ifft(ifftshift(fftshift(fft(Bubble.Contrast,n_pad_fft,2),2).*Greens_function.*phase_shift,2),n_pad_fft,2));
                        %                         pnl.semi_an=pnl.semi_an(:,1:domain.tdimpar);
                        %
                        %                         pnl.semi_anN = sum(pnl.semi_an,1);
                        %                         [RRMS_ErrorG_semi_an(k,m),Rel_ErrorG_semi_an(k,m)] = Error(pnl.semi_anN_cluster,pnl.simN_cluster)   ;
                        
                        %================================================== Compare for passive scattering =======================================
                        %                         passive = strsplit(file.scatterer,'_');
                        %                         if (strcmp(passive{1},'passive'))
                        %                             % The shift due to the time that takes for the planewave to travel to the position of the scatterer
                        %                             % is not taken into account because the output of INCS do not include this
                        %                             % that's why in INCS the movie should be shifted with the z axis, otherwise the result is not logical.
                        %                             % but in case of test, r_from_transducer = sqrt(sum(Bubble.LocGlob(iBubble,3).^2))*1e-3;
                        %                             if (strcmp(passive{2},'lin'))
                        %                                 %  Spectral derivative
                        %                                 R_rigid         =   1E-4;
                        %                                 time_signature  =   (4/3*pi*R_rigid.^3)/(medium.c0^2).*Greens_function;
                        %                                 ddpddt          =   fft(Inc_Press+P_scattered,n_pad_fft,2).*(1i*fftshift(omega,2)).^2 ; % What the scatterer understands
                        %
                        %                             else
                        %  Spectral derivative
                        Amplitude       =   1E-23;
                        time_signature  =   Amplitude.*Greens_function;
                        ddpddt          =   fft((Inc_Press+P_scattered).^2,n_pad_fft,2).*(1i*fftshift(omega,2)).^2 ; % What the scatterer understands
                        
                        %                             end
                        
                        ddpddt = ifftshift(fftshift(ddpddt,2).*time_signature.*phase_shift,2);
                        pnl.an = ddpddt;
                        pnl.an = real(ifft(pnl.an,n_pad_fft,2));
                        if (i>1 && j==1); P_scattered = pnl.an;end % Add the scattered pressure from the other scatterer
                        pnl.anN_cluster = sum(pnl.an(:,1:domain.tdimpar*OSFactor),1);
                        [RRMS_ErrorG_an(k,m),~] = Error(pnl.anN_cluster,pnl.simN_cluster);
                        %
                        %                         end
                        %Safety measure , max 10 plots
                        %                         if ( ((k_end-k_init+1)*(m_end-m_init+1)<10  )  )
                        %                             Plot_Error(tpar_N,pnl,domain,SliceIndex,RRMS_ErrorG_semi_an,RRMS_ErrorG_an,k,m,i,j)
                        %                             set(gca,'FontSize',20)
                        %                             Plot_Spectral_Response(domain,pnl,OSFactor)
                        % %                             text(gcf,'Position', [ 0 0],'String', '(b)','FontSize',35,'FontName','Helvetica')
                        %                             set(gcf,'Color','white')
                        %                             set(gca,'FontSize',20)
                        %                             %% ========================== Save plot ===================================
                        %                             if (strcmp(file.saveplot,'yes'))
                        %                                 %                             saveas(gca,[file.savedir,'/',file.dirname,'_comparison_scattered_pressure_slice',int2str(dslice.num),int2str([i;j;iBubble])','.jpg'])
                        %                                 figure(gcf)
                        %                                 folder = ['../',file.savedir,'/',file.dirname,'_comparison_scattered_pressure_slice',int2str(dslice.num)','_2sc.png'];
                        %                                 cd export_fig/
                        %                                 export_fig  folder -painters -q110
                        %                                 movefile('folder.png',folder)
                        %                                 cd ../
                        %                             end
                        %                         end
                        
                    end
                end
            end
         end
    end
    %% ========================================= PLOT ERROR IN 2D SLICE ==================================================================
    if ( k_end == domain.dimlen(1) && m_end == domain.dimlen(3))
        RRMS_Err_Corr = zeros(domain.dimlen(1),domain.dimlen(3));
        RRMS_Err_Corr(k_init:k_end,m_init:m_end) = RRMS_ErrorG_an(k_init:k_end,m_init:m_end);
        RRMS_Err_Corr(abs(RRMS_Err_Corr-1)<=1E-1)=0;
        RRMS_Err_Corr(RRMS_Err_Corr>=0.2)=0;
        %         RRMS_Err_Corr(RRMS_Err_Corr>0.2)=0;
        figure('units','normalized','Position',[-0.008333333333333,0.028515240904621,0.6953125,0.867256637168142]);
        contourf(domain.par{3},domain.par{1},RRMS_Err_Corr*100)
        BubbleCluster_Colormaps(file)
        hold on
        shading interp;
        xlabel('Z [mm]')
        ylabel('X [mm]')
        c = colorbar;
        c.Label.String = 'RRRMS [%]';
        title(['Analytic vs Simulation, y = ',num2str(dslice.pos(dslice.num)),])
        mean(mean(RRMS_Err_Corr))
        std(std(RRMS_Err_Corr))
        set(gcf,'Color','white')
        set(gca,'FontSize',30)
        if Bubble.N >6
            rectangle('Position',[domain.min_dim(domain.dimval(3)) domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(3))-domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(2))-domain.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor',[190 0 0]/255,'LineWidth',5)
        else
            width = 0.3;
            z1 = domain.par{1}(53);x1=domain.par{3}(75);
            p1 =plot (x1+width*[-1 1 1 -1 -1], z1+width*[-1 -1 1 1 -1],'-','Color',[153 0 0 ]/255,'LineWidth',2) ;
            p1 =plot (x1, z1,'Color',[153 0 0 ]/255) ;
            %         p2 =plot (5.434+width*[-1 1 1 -1 -1], 0.3293+width*[-1 -1 1 1 -1],'-','Color',[204 0 102]/255,'LineWidth',2) ;
            %         p1 =plot (5.434, 0.3293,'Color',[153 0 0 ]/255) ;
            %6.34, -0.494, 58, 77
            %5.434, 0.3293, 68, 66
            plot(Bubble.LocGlob(:,3),Bubble.LocGlob(:,1), 'x','Color',[204 102 0]/255,'MarkerSize',20,'LineWidth',4)
        end
        text('Position', [-1.5 5.5],'String', '(a)','FontSize',35,'FontName','Helvetica')
        %%
        %========================== Save plot ===================================
        if (strcmp(file.saveplot,'yes'))
            figure(gcf)
            %     saveas(f4,[file.savedir,'/',file.dirname,'_3D_cluster_slice',int2str(dslice.num),'.jpg'])
            folder = ['../',file.savedir,'/',file.dirname,'_Comparison_Slice',int2str(dslice.num),'_2sc.png'];
            cd export_fig/
            export_fig  folder -painters -q110
            movefile('folder.png',folder)
            cd ../
        end
    end
    %% ========================================= PLOT Scattered pressure ==================================================================
    
    %     for iBubble = 1:Bubble.N
    %
    %         figure;
    %         ax = gca;
    %         hold on;
    %         plot(domain.tpar,P_scattered(iBubble,:))
    %         ax=gca;
    %         xlabel('time [μsec]')
    %         set(gca,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [μsec]';
    %         ylabel('Pressure [Pa]')
    %         title(['Scattered pressure of Bubble No. ',num2str(Bubbles(iBubble)),' to Bubble No. ' num2str(iBubble), '.'])
    %         hold off;
    %     end
end
end


% Calculate RRMS and Relative Error between analytic and simulation
function [RRMS_ErrorG_an,Rel_ErrorG_an] = Error(p_an,p_sim)

% this is to reduce the very small diffrence in time (dtpar/16)
% max1 = find(max(p_an)==p_an);
% max2 = find(max(p_sim)==p_sim);
% shift = (max2-max1);
% if (shift>0);p_an = [zeros(shift-1,1); p_an(1:end-shift+1)];elseif(shift<0);p_sim = [zeros(-shift-1,1) ;p_sim(1:end+shift+1) ];end
err1 = find(abs(p_sim)>=min(abs(p_sim)));
RRMS_ErrorG_an = sqrt( sum((p_sim(err1) - p_an(err1)).^2,2)./sum(p_an(err1).^2));

% [~,i_maxdiff] = max(abs(p_sim-p_an));
% Rel_ErrorG_an = abs((p_sim(i_maxdiff)-p_an(i_maxdiff))/p_an(i_maxdiff));
Rel_ErrorG_an = 0;
end

%%================================================ Plot Error================================================

function Plot_Error(tpar_N,pnl,domain,SliceIndex,RRMS_ErrorG_semi_an,RRMS_ErrorG_an,k,m,i,j)

f6 = figure('WindowState','maximized');
h_f1 = subplot(2,1,1);
ax = gca; hold on;
p1 = plot(tpar_N ,pnl.simN_cluster,'--o','MarkerIndices',1:50:length(tpar_N),'MarkerSize',4,'color','k','LineWidth',1.5);
% p3 = plot(tpar_N ,pnl.semi_anN_cluster,':x','MarkerIndices',1:50:length(tpar_N),'MarkerSize',4,'color',[17 17 17]./255);
p2 = plot(tpar_N ,pnl.anN_cluster,'-d','MarkerIndices',1:50:length(tpar_N),'MarkerSize',4,'color',[17 17 17]./255,'LineWidth',1.5);
p3= plot([],[]);
h =[p1;p2;p3];

ax.YLabel.String = 'Pressure [Pa]';
legendBubble{2} =['Simulation, RRMSE =',num2str(round(RRMS_ErrorG_an(k,m) *1000)/10),'[%]'];
legendBubble{1} =['Analytic '];
% legendBubble{3} =['Semi Analytic'];
% xlim([min(tpar_N) max(tpar_N)])
% RRMS_semi_an = ['RRMS (Semi) ='        ,num2str(round(RRMS_ErrorG_semi_an(k,m) *1000)/10),'[%]'];
% RelE_semi_an = ['|Rel.Error| (Semi) = ',num2str(round(Rel_ErrorG_semi_an*1000)/10),'[%]'];
RRMS_an      = ['RRMS (Full) ='        ,num2str(round(RRMS_ErrorG_an(k,m) *1000)/10),'[%]'];
% RelE_an      = ['|Rel.Error| (Full) = ',num2str(round(Rel_ErrorG_an*1000)/10),'[%]'];
% annotation('textbox', [0.75, 0.81, 0.1, 0.03], 'String', [RRMS_an newline RRMS_semi_an])
xlim([4e-6 10e-6])
legend(h,legendBubble{:},'Location','northwest')
title(['Analytic vs Simulation of MB scattered pressure at (x,y,z) = (',...
    num2str(round(100*domain.par{domain.dimval(2)}(SliceIndex(1)))/100),', ',...
    num2str(round(100*domain.par{domain.dimval(1)}(SliceIndex(2)))/100),', '...
    num2str(round(100*domain.par{domain.dimval(3)}(SliceIndex(3)))/100),') [mm]'])

set(ax,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [microsec]';
end

%%==================================== Plot Spectral Response ====================================================

function Plot_Spectral_Response(domain,pnl,OSFactor)

F_s = 1/(domain.dtpar/OSFactor*domain.PPW_t);                          % [Hz] , 1/dt
[f_p,Sp_psim] = Freq_Calc(pnl.simN_cluster,F_s);            % Spectrum of simulation
[~,Sp_pan] = Freq_Calc(pnl.anN_cluster,F_s);                % Spectrum of analytical pulse
[~,Sp_psemi_an] = Freq_Calc(pnl.semi_anN_cluster,F_s);      % Spectrum of semi_analytical pulse

% Create plot of frequency spectrum for driving pulse and radial oscillations
%                     f7 = figure('WindowState','maximized');
h_f2 = subplot(2,1,2);
ax = gca;   hold on;
m_plot=plot(f_p(1:domain.tdimpar),Sp_psim(1:domain.tdimpar)-max(Sp_psim), 'k--', 'LineWidth', 1.5);
fig_color = get(m_plot,'Color');
pr=plot(f_p(1:domain.tdimpar),Sp_pan(1:domain.tdimpar)-max(Sp_pan),'k-','LineWidth', 1.5);
% psr=plot(f_p(1:domain.tdimpar),Sp_psemi_an(1:domain.tdimpar)-max(Sp_psemi_an),':','LineWidth', 1, 'color', fig_color);
legend([m_plot, pr], 'Simulation', 'Analytic')
ylim([-45 max(Sp_pan(1:domain.tdimpar)-max(Sp_pan))]);
ylabel('normalized amplitude [dB]');
% title('Frequency Spectrum');
box on;
xlim([0 F_s/OSFactor]);

set(gca,'XTickLabel',ax.XTick*1e-6) ; ax.XLabel.String = 'Frequency [MHz]';
hold off;
end

% Calculate the frequency spectrum based on a sampling frequency and the pressure values
function [f_p,S_p] = Freq_Calc(Pres,F_s)

Nfft_p = 2^nextpow2(2*length(Pres));                      % Zero-padding to improve fft performance until length equal to 2N
Y_p = abs(fft(Pres,Nfft_p));                              % Convert the driving pulse to the frequency domain
f_p = linspace(0,1,Nfft_p/2+1)*F_s;                     % [Hz] Frequency sample Number
% Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
% Other way to define frequency domain : f_p = Fs*(0:(Nfft_p/2))/Nfft_p
Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
S_p = 20*log10(Y_p(Np));                                     % Convert pressure to amplitude [db]

end
