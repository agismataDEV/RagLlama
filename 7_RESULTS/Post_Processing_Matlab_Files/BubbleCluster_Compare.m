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


function BubbleCluster_Compare_Individually(medium,domain,dslice,file,plin,pnl,Bubble)
disp(' ')
disp('Calculating the difference between analytic and simulation solution...')
if (Bubble.N <=10 && Bubble.N>0)
    
    Sim_OSFactor=1;
    %%
    i_start = 1;
    i_end = 200;
    depth =  domain.dimlen(3)-20;
    figure('WindowState','maximized');
    hold on;  plot(domain.tpar(i_start:i_end)*1E6,plin.data(i_start:i_end,round(size(plin.data,3)/2)+1,depth)*1E-3,'-','LineWidth',2)
    hold on;  plot(domain.tpar(i_start:i_end)*1E6,pnl.anN_cluster(i_start:i_end,round(size(plin.data,3)/2)+1,depth)*1E-3, '--','LineWidth',2)
    xlabel('Time [usec]')
    ylabel('Pressure [kPa]')
    set(gca,'FontSize',20)
    set(gcf,'Color','white')
    legend('Linear','With MBs')
    grid on;
    title(['INCS Results after 10 Iterations @ (X,Y,Z) = (', num2str(domain.par{1}(70)),', ',num2str(0),', ', num2str(domain.par{3}(depth)),') [mm]'])
    
    F_s = 1/(2*domain.dtpar);                          % [Hz] , 1/dt
    [f_p,Sp_p] = Freq_Calc(plin.data(i_start:i_end,70,depth),1/(2*domain.dtpar));            % Spectrum of simulation
    figure('WindowState','maximized');
    hold on;  plot(f_p*1E-6,Sp_p,'-','LineWidth',2)
    [f_s,Sp_s] = Freq_Calc(pnl.contrastdata(i_start:i_end,70,depth),1/(2*domain.dtpar));            % Spectrum of simulation
    hold on;  plot(f_s*1E-6,Sp_s,'--','LineWidth',2)
    xlabel('Frequency [MHz]')
    ylabel('Amplitude [dB]')
    set(gca,'FontSize',20)
    set(gcf,'Color','white')
    legend('Linear','With MBs')
    grid on;
%     ylim([max([Sp_p; Sp_s])-70 max([Sp_p; Sp_s])])
    %% ============== Comparison of time signature based on ODE Solvers ====================================
    if (strcmp(file.scatterer,'active'))
        %%
        % Load pressure values from file
        %             filename_pres = filename_pres(1:end-1);
        filename_pres = split(cellstr(ls([file.dirname '/Bubbles/' 'v_dd_norm_004*.txt'])));
        for ifile =1:1
            fileID_pres = fopen([file.dirname '/Bubbles/' filename_pres{ifile} ],'r');
            size_pres = [1 Inf];
            formatSpec = '%f';
            
            V_dd_sim_loaded = fscanf(fileID_pres,formatSpec,size_pres);
            
            n_pad = (length(V_dd_sim_loaded)/domain.tdimpar-2)/2;
            n_pad = 4;
            iBubble_ofinterest =  1 ;
            if (length(V_dd_sim_loaded)>(1+2*domain.tdimpar+2*n_pad*domain.tdimpar)); iBubble_ofinterest = 1;end
            V_dd_sim = V_dd_sim_loaded((iBubble_ofinterest-1)*(2*domain.tdimpar+2*n_pad*domain.tdimpar) + 1 :iBubble_ofinterest*(2*domain.tdimpar+2*n_pad*domain.tdimpar));
            
            V_dd_Time = V_dd_sim(1:domain.tdimpar);
            V_dd_Pres = V_dd_sim(domain.tdimpar + 1 : 2*domain.tdimpar) * medium.rho0;
            V_dd_TimePad = V_dd_sim(2*domain.tdimpar + 1:(n_pad+2)*domain.tdimpar);
            V_dd_PresPad = V_dd_sim((n_pad+2)*domain.tdimpar + 1 : (n_pad+2)*domain.tdimpar+n_pad*domain.tdimpar)* medium.rho0;
            
            figure('WindowState','maximized');
            hold on; plot(V_dd_Time*1E6,V_dd_Pres,V_dd_TimePad*1E6,V_dd_PresPad)
            legend('Initial Pulse', 'Resampled Pulse')
            set(gca,'FontSize',20)
            set(gcf,'Color','white')
            xlabel('Time [usec]')
            ylabel('Contrast Source term [Pa/m^2]')
            
            F_s = 1/(2*domain.dtpar/n_pad);                          % [Hz] , 1/dt
            [f_p,Sp_p] = Freq_Calc(V_dd_Pres,1/(2*domain.dtpar));            % Spectrum of simulation
            [f_ppad,Sp_ppad] = Freq_Calc(V_dd_PresPad,F_s);                % Spectrum of analytical pulse
            figure('WindowState','maximized');
            hold on;  plot(f_p(1:domain.tdimpar)*1E-6,Sp_p(1:domain.tdimpar)-max(Sp_p),f_ppad(1:domain.tdimpar)*1E-6,Sp_ppad(1:domain.tdimpar)-max(Sp_ppad))
            legend('Initial Pulse', 'Resampled Pulse')
            set(gca,'FontSize',20)
            set(gcf,'Color','white')
            xlabel('Frequency [MHz]')
            ylabel('Normalized amplitude [dB]')
            
            % Create Simulation Contrast, correct if n_padded in INCS code
            BubbleTimeSim = resampleSINC(domain.tpar,n_pad);
            Sim_OSFactor = n_pad*length(Bubble.Contrast(1,:))/length(BubbleTimeSim);
            if Sim_OSFactor>1 ; BubbleCon_Sim = downsample(Bubble.Contrast(1,:),Sim_OSFactor/n_pad); else ; BubbleCon_Sim = resampleSINC(Bubble.Contrast(1,:),n_pad/Sim_OSFactor); end
        end
        fclose('all')
        %%
        filename_pres = split(cellstr(ls([file.dirname '/Bubbles/output_file_pressure_001*.txt'])));
        for ifile = 1:1
            fileID_pres = fopen([file.dirname '/Bubbles/' filename_pres{ifile} ],'r');
            size_pres = [1 Inf];
            formatSpec = '%f';
            
            BubbleVal_loaded = fscanf(fileID_pres,formatSpec,size_pres);
            
            n_pad = ((length(BubbleVal_loaded)-1)/domain.tdimpar-2)/2;
            n_pad = 4;
            iBubble_ofinterest = 1;
            if (length(BubbleVal_loaded)>(1+2*domain.tdimpar+2*n_pad*domain.tdimpar)); iBubble_ofinterest = 1;end
            BubbleVal = BubbleVal_loaded((iBubble_ofinterest-1)*(1+2*domain.tdimpar+2*n_pad*domain.tdimpar) + 2 :iBubble_ofinterest*(1+2*domain.tdimpar+2*n_pad*domain.tdimpar));
            
            BubbleTime = BubbleVal(1:domain.tdimpar);
            BubblePres = BubbleVal(domain.tdimpar + 1 : 2*domain.tdimpar);
            BubbleTimePad = BubbleVal(2*domain.tdimpar + 1:(n_pad+2)*domain.tdimpar);
            BubblePresPad = BubbleVal((n_pad+2)*domain.tdimpar + 1 : (n_pad+2)*domain.tdimpar+n_pad*domain.tdimpar);
            figure('WindowState','maximized');
            hold on; plot(BubbleTime*1E6,BubblePres,BubbleTimePad*1E6,BubblePresPad)
            legend('Initial Pulse', 'Resampled Pulse')
            set(gca,'FontSize',20)
            set(gcf,'Color','white')
            xlabel('Time [usec]')
            ylabel('Pressure [Pa]')
            grid on;
            
            F_s = 1/(2*domain.dtpar/n_pad);                          % [Hz] , 1/dt
            [f_p,Sp_p] = Freq_Calc(BubblePres,1/(2*domain.dtpar));            % Spectrum of simulation
            [f_ppad,Sp_ppad] = Freq_Calc(BubblePresPad,F_s);                % Spectrum of analytical pulse
            figure('WindowState','maximized');
            hold on;  plot(f_p(1:domain.tdimpar)*1E-6,Sp_p(1:domain.tdimpar)-max(Sp_p),f_ppad(1:domain.tdimpar)*1E-6,Sp_ppad(1:domain.tdimpar)-max(Sp_ppad))
            legend('Initial Pulse', 'Resampled Pulse')
            set(gca,'FontSize',20)
            set(gcf,'Color','white')
            xlabel('Frequency [MHz]')
            ylabel('Normalized amplitude [dB]')
            grid on;
            
            % Create Simulation Contrast, correct if n_padded in INCS code
            BubbleTimeSim = resampleSINC(domain.tpar,n_pad);
            Sim_OSFactor = n_pad*length(Bubble.Contrast(1,:))/length(BubbleTimeSim);
            if Sim_OSFactor>1 ; BubbleCon_Sim = downsample(Bubble.Contrast(1,:),Sim_OSFactor/n_pad); else ; BubbleCon_Sim = resampleSINC(Bubble.Contrast(1,:),n_pad/Sim_OSFactor); end
        end
        fclose('all')
        %%
%         BubblePresPad = resampleSINC(plin.data(:,72,80)',n_pad);
        BubbleTimeAn = BubbleTimePad;
        V_dd_norm = Bubble_SimAnalytic(BubblePresPad,BubbleTimeAn,medium.c0,medium.rho0,medium.freq0);         % Solution of RP Solver with ode45 MATLAB solver
        BubbleCon_an = V_dd_norm*medium.rho0;                                                             % Contrast Source Term
        %========================== Save plot ===================================
        if (strcmp(file.saveplot,'yes'))
            saveas(gcf,[file.savedir,'/',file.dirname,'_bubble_motion_slice',int2str(dslice.num),'.jpg'])
        end
        %
        f5 = figure('WindowState','maximized');
        ax = gca; hold on;
        plot(BubbleTimeAn*1E6,V_dd_PresPad,'--','color','k');
        plot(BubbleTimeAn*1E6,BubbleCon_an,'-','color',[17 17 17]./255);
        legend('Simulation', 'Matlab ODE')
        set(gca,'FontSize',20)
        set(gcf,'Color','white')
        xlabel('Time [usec]')
        ylabel('Contrast Source term [Pa/m^2]')
        grid on;
        
        F_s = 1/(domain.dtpar/n_pad);                          % [Hz] , 1/dt
        [f_p,Sp_p] = Freq_Calc(BubbleCon_an,F_s);            % Spectrum of simulation
        [~,Sp_ppad] = Freq_Calc(V_dd_PresPad,F_s);                % Spectrum of analytical pulse
        f_p = f_p*1E-6/2;
        figure('WindowState','maximized');plot(f_p(1:domain.tdimpar),Sp_p(1:domain.tdimpar)-max(Sp_p));hold on;plot(f_p(1:domain.tdimpar),Sp_ppad(1:domain.tdimpar)-max(Sp_ppad))
        set(gca,'FontSize',20)
        set(gcf,'Color','white')
        xlabel('Frequency [MHz]')
        ylabel('Normalized amplitude [dB]')
        legend('Simulation', 'Matlab ODE')
        grid on;
        fclose('all')
    end
    
    %% =============================================== Compare simulation with analutical Results =======================================================
    % This is working for 1 or 2 Bubbles and up to 1st order multiple scattering
    % it can easily expanded to more but these are enough for the validation of the code
    %         close all
    
    OSFactor = 	1;
    R_rigid = repmat(5E-6,Bubble.N,1);
    time_signature  =   (4/3*pi*R_rigid.^3)/(medium.c0^2);
    medium.rho1     =   10;
    medium.c1       =   100;
    time_signature  =   (4/3*pi*R_rigid.^3)*medium.rho0.*(1./(medium.rho1*medium.c1.^2)- 1/(medium.rho0*medium.c0^2));
    
    %======== Driving Frequency ======================
    freq        = medium.freq0;
    w_driving   = 2*pi*freq;
    gauss_dl    = 3/freq;
    gauss_win   = 1.5/freq;
    
    domain.tpar_OS = interp(domain.tpar,OSFactor);
    domain.dtpar_OS = domain.tpar_OS(2) - domain.tpar_OS(1);
    domain.tdimpar_OS = length(domain.tpar_OS);
    
    Inc_Press   = repmat(1E5*exp( - ((domain.tpar_OS-gauss_dl)/(gauss_win/2)).^2).* sin(w_driving*(domain.tpar_OS-gauss_dl)),Bubble.N,1);
    Fs          = 1/domain.dtpar_OS;  dF = Fs/domain.tdimpar_OS;
    if mod(domain.tdimpar_OS,2) == 0 % Nt is even
        f = -Fs/2 : dF : Fs/2-dF;
    else            % Nt is odd
        f = [sort(-1*(dF:dF:Fs/2)) (0:dF:Fs/2)];
    end
    omega       = 2*pi*f;
    omega       = [omega(omega>=0) omega(omega<0)];
    omega       = repmat(omega,Bubble.N,1);
    
    % This is only for the case of planewave. This is a time shift for the travelling of a planewave.
    Inc_Press = real(ifft(fft(Inc_Press,[],2).*exp(-1i*omega/medium.c0.*Bubble.LocGlob(:,3)*1E-3),[],2));
    % Make denser the matrix to increase accuracy
    pdiff  = pnl.contrastdata;  % If you do not save contrastdata,  pdiff = plin.data - pnl.data
    
    %======= Initialize Arrays ========================
    pnl.simN                = zeros(domain.tdimpar*OSFactor,1);
    pnl.anN                 = zeros(domain.tdimpar*OSFactor,1);
    pnl.semi_anN            = zeros(domain.tdimpar*OSFactor,1);
    pnl.anN_map             = zeros(domain.tdimpar*OSFactor,1);
    
    [~,SliceIndex(2)] = min(abs(dslice.pos(dslice.num)-domain.par{domain.dimval(1)}));
    % Iterate over k in x axis and m in z axis. This implementation is done for a y slice so the index in y equals always to the closest point of the slice position (SliceIndex(2))
    %% Generate the pressure that each microbubble understands at the iteration of interest
    P_scattered = zeros(Bubble.N,domain.tdimpar_OS);
    TotalPress = Inc_Press;
    for j = 1: domain.beamiterations-1
        P_scattered = zeros(Bubble.N,domain.tdimpar_OS);
        for iBubble = 1:ceil(Bubble.N/2)
            r_from_scatterer = sqrt((Bubble.LocGlob(iBubble,1) - Bubble.LocGlob(:,1)).^2 + (Bubble.LocGlob(iBubble,2) - Bubble.LocGlob(:,2)).^2 + (Bubble.LocGlob(iBubble,3) - Bubble.LocGlob(:,3)).^2)*1e-3;
            Greens_function = (exp(-1i*omega.*r_from_scatterer/medium.c0)./(4*pi*r_from_scatterer)); % 1/(4πr) and the time shift due to wave travel
            Greens_function(isinf(Greens_function)) = 0;
            
            ddpddt          =   repmat(fft(TotalPress(iBubble,:),[],2).*(1i*omega(iBubble,:)).^2,Bubble.N,1) ; % What the scatterer understands
            ddpddt = ddpddt.* time_signature(iBubble) .* Greens_function;
            pnl.an = ddpddt;
            pnl.an = real(ifft(pnl.an,[],2));
            
            iBubble_end = Bubble.N - iBubble+1;
            r_from_scatterer_end = sqrt((Bubble.LocGlob(iBubble_end,1) - Bubble.LocGlob(:,1)).^2 + (Bubble.LocGlob(iBubble_end,2) - Bubble.LocGlob(:,2)).^2 + (Bubble.LocGlob(iBubble_end,3) - Bubble.LocGlob(:,3)).^2)*1e-3;
            Greens_function_end = (exp(-1i*omega.*r_from_scatterer_end/medium.c0)./(4*pi*r_from_scatterer_end)); % 1/(4πr) and the time shift due to wave travel
            Greens_function_end(isinf(Greens_function_end)) = 0;
            
            ddpddt_end          =   repmat(fft(TotalPress(iBubble_end,:),[],2).*(1i*omega(iBubble_end,:)).^2,Bubble.N,1) ; % What the scatterer understands
            ddpddt_end = ddpddt_end.* time_signature(iBubble_end) .* Greens_function_end;
            pnl.an_end = ddpddt_end;
            pnl.an_end = real(ifft(pnl.an_end,[],2));
            
            % This array contains for each row, the pressure that this bubble understands as a summation of all its neighbors
            % So the first row is for
            P_scattered = P_scattered + pnl.an + pnl.an_end;
        end
        TotalPress      =   Inc_Press + P_scattered;
        
    end
    
    %%
    k_init      =   1;%79;
    k_end       =   domain.dimlen(1);
    m_init      =   1;%domain.dimlen(3);
    m_end       =   domain.dimlen(3);
    % % %
%     %
%     m_init = idx1;
%     k_init = idx3;
    m_init      =   69;49;77;
    m_end       =   m_init;
    k_init      =   10;11;58;
    k_end       =   k_init;
    %     %
    RRMS_ErrorG_semi_an = zeros(k_end,m_end);
    RRMS_ErrorG_an = zeros(k_end,m_end);
    for i = domain.beamiterations:domain.beamiterations
        TotalPress_lastiter      =   TotalPress;
        ddpddt          =   fft(TotalPress_lastiter,[],2).*(1i*omega).^2 ; % What the scatterer understands
%         syms t
%         p_expr  = 1E4*exp( - ((t-gauss_dl)/(gauss_win/2)).^2).* sin(w_driving*(t-gauss_dl));
%         p = diff(p_expr,t,2);p=strrep(char(p),'*','.*');p=strrep(char(p),'^','.^');p=strrep(char(p),'/','./');p =strrep(p,'t','domain.tpar_OS');
%         ddpddt =fft((((800000.*exp(-(4.*(gauss_dl - domain.tpar_OS).^2)./gauss_win.^2).*sin(w_driving.*(gauss_dl - domain.tpar_OS)))./gauss_win.^2 + 100000.*w_driving.^2.*exp(-(4.*(gauss_dl - domain.tpar_OS).^2)./gauss_win.^2).*sin(w_driving.*(gauss_dl - domain.tpar_OS)) - (1600000.*exp(-(4.*(gauss_dl - domain.tpar_OS).^2)./gauss_win.^2).*sin(w_driving.*(gauss_dl - domain.tpar_OS)).*(2.*gauss_dl - 2.*domain.tpar_OS).^2)./gauss_win.^4 + (800000.*w_driving.*exp(-(4.*(gauss_dl - domain.tpar_OS).^2)./gauss_win.^2).*cos(w_driving.*(gauss_dl - domain.tpar_OS)).*(2.*gauss_dl - 2.*domain.tpar_OS))./gauss_win.^2)));
%         ddpddt = ddpddt.*exp(-1i*omega/medium.c0.*Bubble.LocGlob(:,3)*1E-3);
         %======= Initialize Arrays ========================
        % Initialize for higher order multiple scattering . This is done because on previous j's scattered pressure is calculated
        pnl.simN_cluster        = zeros(domain.tdimpar*OSFactor,1);
        pnl.anN_cluster         = zeros(domain.tdimpar*OSFactor,domain.dimlen(1),domain.dimlen(3));
        pnl.semi_anN_cluster    = zeros(domain.tdimpar*OSFactor,domain.dimlen(1),domain.dimlen(3));
        for k = k_init: k_end
            for m = m_init:m_end
                % Iterate to include higher orders of multiple scattering
                % Find the position of interest which is the position of slice
                % Find the grid points with which the bubble has the smallest distance
                SliceIndex(1) = k ;  SliceIndex(3) = m;
                %================================================== GLOBAL VARIABLES =============================================================
                % This is used to show the propagation in the horizontal axis,This is because of the results of INCS
                % As in the movies, there should be a shift in
                r_from_scatterer = sqrt((Bubble.LocGlob(:,domain.dimval(2)) - domain.par{domain.dimval(2)}(SliceIndex(1))).^2+(Bubble.LocGlob(:,domain.dimval(1)) - dslice.pos(dslice.num)).^2 ...
                    +(Bubble.LocGlob(:,domain.dimval(3))-domain.par{domain.dimval(3)}(SliceIndex(3))).^2)*1e-3;
                z_shift = (Bubble.LocGlob(:,domain.dimval(3))-domain.par{domain.dimval(3)}(SliceIndex(3)) )*1e-3;
                
                phase_shift = 1;
                Greens_function = (exp(-1i*omega.*r_from_scatterer/medium.c0)./(4*pi*r_from_scatterer)); % 1/(4πr) and the time shift due to wave travel
                Greens_function(isinf(Greens_function)) = 0;
                
                %================================================== Compare for passive scattering =======================================
                passive = strsplit(file.scatterer,'_');
                %                             if (strcmp(passive{1},'passive'))
                % The shift due to the time that takes for the planewave to travel to the position of the scatterer
                % is not taken into account because the output of INCS do not include this
                % that's why in INCS the movie should be shifted with the z axis, otherwise the result is not logical.
                % but in case of test, r_from_transducer = sqrt(sum(Bubble.LocGlob(iBubble,3).^2))*1e-3;
                
                %                                 if (strcmp(passive{2},'lin'))
                %Spectral derivative
             
                ddpddt_lastiter = ddpddt.* time_signature .* Greens_function;
                pnl.an = ddpddt_lastiter;
                pnl.an = real(ifft(pnl.an,[],2));
                pnl.anN_cluster(:,k,m) = sum(pnl.an,1)';
                
                %                             elseif (strcmp(passive{1},'active'))
                %                                 ddpddt = ifftshift(fftshift(fft(Bubble.Contrast(iBubble,:))) .* Greens_function);
                %                                 pnl.an = ddpddt;
                %                                 pnl.anN = real(ifft(ddpddt))';%interp(pnl.an,OSFactor)';
                %
                %                                 pnl.anN_cluster = pnl.anN_cluster + pnl.anN;
                %                             end
                %                             pnl.anN_map(:,k,m) = pnl.anN_cluster;
                
                %================================================== SEMI-ANALYTIC RESULTS SHIFTED =============================================================
                % This phase shift is due to the shift that was implemented in the incident pressure field in the beginning
                % When a comoving time window is used then all the results need to be shifted in the correct position
                if (domain.a_t); phase_shift = exp(-1i*omega.*repmat(Bubble.LocGlob(:,3)*1E-3,1,domain.tdimpar_OS)/medium.c0); end
                pnl.semi_an = Bubble.Contrast ; % Semi - analytic because checking only the Green's function
                pnl.semi_an = real(ifft(fft(pnl.semi_an,[],2).*Greens_function.*phase_shift,[],2));
                pnl.semi_anN_cluster(:,k,m) = sum(pnl.semi_an,1)';%interp(pnl.semi_an,OSFactor)'; % Normalize to calculate error
                
                %================================================== SIMULATION RESULTS =============================================================
                pnl.sim = squeeze(pdiff(:,SliceIndex(1),SliceIndex(3)));
                if (domain.a_t);pnl.sim = real(ifft(fft(pnl.sim').*phase_shift(1,:).* exp(1i*omega(1,:).*z_shift(1,1)/medium.c0)))';end
                % The check for the time signature for microbubbles is done from the previous step of the comparison of the ode solvers
                pnl.simN_cluster = pnl.sim;%interp(pnl.sim,OSFactor); % Normalize to calculate error
                
                [RRMS_ErrorG_an(k,m),Rel_ErrorG_an(k,m)] = Error(pnl.anN_cluster(:,k,m),pnl.simN_cluster);
                [RRMS_ErrorG_semi_an(k,m),Rel_ErrorG_semi_an(k,m)] = Error(pnl.semi_anN_cluster(:,k,m),pnl.simN_cluster);
                
                %Safety measure , max 10 plots
                if ( (k_end-k_init+1)*(m_end-m_init+1)<10)
                    Plot_Error(domain.tpar_OS,pnl,domain,SliceIndex,RRMS_ErrorG_semi_an,RRMS_ErrorG_an,k,m,Bubble.N,i)
                    Plot_Spectral_Response(domain,pnl,OSFactor,k,m)
                    %% ========================== Save plot ===================================
                    if (strcmp(file.saveplot,'yes'))
                        saveas(gca,[file.savedir,'/',file.dirname,'_comparison_scattered_pressure_slice',int2str(dslice.num),int2str([i;j;iBubble])','.jpg'])
                        %saveas(f7,[file.savedir,'/',file.dirname,'_frequencyspectrum',int2str(dslice.num),int2str([i;j])','.jpg'])
                    end
                end
            end
        end
    end
    
    %% ========================================= PLOT ERROR IN 2D SLICE ==================================================================
    %     if ( k_end == domain.dimlen(1) && m_end == domain.dimlen(3))
    RRMS_Err_Corr = RRMS_ErrorG_an;
    RRMS_Err_Corr(abs(RRMS_Err_Corr)>0.9)=0;
%         RRMS_Err_Corr(abs(RRMS_Err_Corr)>0.2)=0;
    figure('WindowState','maximized');
    ax1 = axes;
    length_z = 1;
    length_x = 1;
    contourf(domain.par{3}(length_z:end-length_z-1),domain.par{1}(length_x:end-length_x-1),RRMS_Err_Corr(length_x:end-length_x-1,length_z:end-length_z-1)*100);
%     contourf(RRMS_Err_Corr*100)
    xlabel('Z [mm]')
    ylabel('X [mm]')
    c1 = colorbar(ax1);
    c1.Label.String = 'RRRMS [%], Estimator : Peak Pressue';
    set(gca,'FontSize',20)
    
    %         ax2 = axes;
    %         c2 = colorbar(ax2);
    %         c2.Label.String = 'RRRMS [%]';
    %         axis off;
    %     if (Bubble.N >6)
    %         %         rectangle('Position',[domain.min_dim(domain.dimval(3)) domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(3))-domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(2))-domain.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor',[190 0 0]/255,'LineWidth',5)
    %         rectangle('Position',[domain.min_dim(domain.dimval(3)) domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(3))-domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(2))-domain.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor','white','LineWidth',5)
    %     else
    %         plot(ax1,Bubble.LocGlob(:,3),Bubble.LocGlob(:,1), 'x','Color','white','MarkerSize',20,'LineWidth',4)
    %     end
    BubbleCluster_Colormaps(file)
    shading interp;
    
    title(['Analytic vs Simulation of MB scattered pressure at y = ',num2str(dslice.pos(dslice.num)),])
    mean(mean(RRMS_Err_Corr))
    std(std(RRMS_Err_Corr))
    set(gca,'FontSize',20)
    set(gcf,'Color','white');
    %% ========================================= PLOT Scattered pressure ==================================================================
    
    for iBubble = 1:Bubble.N
        
        figure;
        ax = gca;
        hold on;
        plot(domain.tpar_OS,P_scattered(iBubble,:))
        ax=gca;
        xlabel('time [μsec]')
        set(gca,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [μsec]';
        ylabel('Pressure [Pa]')
        title(['Scattered pressure of Bubble No. ',num2str(Bubbles(iBubble)),' to Bubble No. ' num2str(iBubble), '.'])
        hold off;
    end
end
end

%% Calculate RRMS and Relative Error between analytic and simulation
function [RRMS_ErrorG_an,Rel_ErrorG_an] = Error(p_an,p_sim)

% this is to reduce the very small diffrence in time (dtpar/16)
% max1 = find(max(p_an)==p_an);
% max2 = find(max(p_sim)==p_sim);
% shift = (max2-max1);
% if (shift>0);p_an = [zeros(shift-1,1); p_an(1:end-shift+1)];elseif(shift<0);p_sim = [zeros(-shift-1,1) ;p_sim(1:end+shift+1) ];end
err1 = find(abs(p_sim)>=min(abs(p_sim)));
meanval = 0;max(p_an(err1));
RRMS_ErrorG_an = sqrt( sum((p_sim(err1) - p_an(err1)).^2,1)./sum((p_an(err1)-meanval).^2 ));

[~,i_maxdiff] = max(abs(p_sim-p_an));
% Rel_ErrorG_an = abs((p_sim(i_maxdiff)-p_an(i_maxdiff))/p_an(i_maxdiff));
% Rel_ErrorG_an = sum((p_srim(err1) - p_an(err1))/length(err1)); % MAE
% Rel_ErrorG_an = sqrt(sum((p_sim(err1) - p_an(err1)).^2/length(err1)));
Rel_ErrorG_an = sqrt( sum((p_sim(err1) - p_an(err1)).^2,1)./sum((p_an(err1)-max(p_an(err1))).^2 )); % Max as an estimator


end

%% ================================================ Plot Error================================================

function Plot_Error(tpar_N,pnl,domain,SliceIndex,RRMS_ErrorG_semi_an,RRMS_ErrorG_an,k,m,iBubble,i)
f6 = figure('WindowState','maximized');
h_f1 = subplot(2,1,1);
ax = gca; hold on;
p1 = plot(tpar_N*1E6 ,pnl.simN_cluster,'--o','MarkerIndices',1:50:length(tpar_N),'MarkerSize',4,'color','k');
p3 = plot(tpar_N*1E6 ,pnl.semi_anN_cluster(:,k,m),':x','MarkerIndices',1:50:length(tpar_N),'MarkerSize',4,'color',[17 17 17]./255);
p2 = plot(tpar_N*1E6 ,pnl.anN_cluster(:,k,m),'-d','MarkerIndices',1:50:length(tpar_N),'MarkerSize',4,'color',[17 17 17]./255);
h =[p1;p2;p3];

legendBubble{iBubble,1} =['Simulation'];
legendBubble{iBubble,2} =['Analytic, Bubble No.',int2str(iBubble)];
legendBubble{iBubble,3} =['Semi Analytic, Bubble No.',int2str(iBubble)];
legend(h,legendBubble{iBubble,:},'Location','northwest')
title(['Analytic vs Simulation of MB scattered pressure at (x,y,z) = (',...
    num2str(round(100*domain.par{domain.dimval(2)}(SliceIndex(1)))/100),', ',...
    num2str(round(100*domain.par{domain.dimval(1)}(SliceIndex(2)))/100),', '...
    num2str(round(100*domain.par{domain.dimval(3)}(SliceIndex(3)))/100),') [mm] , iteration ',num2str(i),', contribution of ',num2str(iBubble),' Scatterer(s)'])
xlim([tpar_N(1) tpar_N(end)]*1E6)
RRMS_semi_an = ['RRMS (Semi) ='        ,num2str(round(RRMS_ErrorG_semi_an(k,m) *1000)/10),'[%]'];
% RelE_semi_an = ['|Rel.Error| (Semi) = ',num2str(round(Rel_ErrorG_semi_an*1000)/10),'[%]'];
RRMS_an      = ['RRMS (Full) ='        ,num2str(round(RRMS_ErrorG_an(k,m) *1000)/10),'[%]'];
% RelE_an      = ['|Rel.Error| (Full) = ',num2str(round(Rel_ErrorG_an*1000)/10),'[%]'];
annotation('textbox', [0.75, 0.81, 0.1, 0.1], 'String', [RRMS_an newline RRMS_semi_an])
xlabel('Time [usec]')
ylabel('Pressure [Pa]')
set(gca,'FontSize',20)
set(gcf,'color','white')
grid on;

end
%%
%%==================================== Plot Spectral Response ====================================================

function Plot_Spectral_Response(domain,pnl,OSFactor,k,m)

F_s = 1/(domain.dtpar/OSFactor*2);                          % [Hz] , 1/dt
[f_p,Sp_psim] = Freq_Calc(pnl.simN_cluster,F_s);            % Spectrum of simulation
[~,Sp_pan] = Freq_Calc(pnl.anN_cluster(:,k,m),F_s);                % Spectrum of analytical pulse
[~,Sp_psemi_an] = Freq_Calc(pnl.semi_anN_cluster(:,k,m),F_s);      % Spectrum of semi_analytical pulse

% Create plot of frequency spectrum for driving pulse and radial oscillations
%                     f7 = figure('WindowState','maximized');
h_f2 = subplot(2,1,2);
ax = gca;   hold on;
m_plot=plot(f_p(1:domain.tdimpar_OS)*1E-6,Sp_psim(1:domain.tdimpar_OS)-max(Sp_psim), '--', 'LineWidth', 1);
fig_color = get(m_plot,'Color');
pr=plot(f_p(1:domain.tdimpar_OS)*1E-6,Sp_pan(1:domain.tdimpar_OS)-max(Sp_pan),'-','LineWidth', 1, 'color', fig_color);
psr=plot(f_p(1:domain.tdimpar_OS)*1E-6,Sp_psemi_an(1:domain.tdimpar_OS)-max(Sp_psemi_an),':','LineWidth', 1, 'color', fig_color);
legend([m_plot, pr,psr], 'Simulation', 'Analytic Pulse', 'Semi-Analytic Pulse' )
ylim([-85 max(Sp_pan(1:domain.tdimpar_OS)-max(Sp_pan))]);
ylabel('normalized amplitude [dB]'); title('power spectra'); box on;
xlabel('Frequency [MHz]')
xlim([0 F_s*1E-6/OSFactor]);

hold off;
set(gca,'FontSize',20)
set(gcf,'color','white')
grid on;

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