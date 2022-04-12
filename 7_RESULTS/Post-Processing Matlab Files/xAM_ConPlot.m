% CREATING MATRICES FOR NONLINEAR AND CONTRAST DATA AND PLOTTING THE PRESSURE FIELD
%
% [pnl] = xAM_ConPlot(medium,domain,slice,file,plin,pnl,Bubble)
% Create the array for the linear pressure field using a butterworth filter of 8th order.
% Create the arrays for all the harmonics up to the sampling frequency used for the discretization
% In case of 1 scatterer due to the small pressure field , a hard bandpass filter is used
% in the frequency domain and the blackmann window in the time domain.
%
% Inputs:
% medium        struct      Contains all the medium properties, e.g Density, Speed of Sound etc
% domain        struct      Contains all the computational domain properties, e.g space and time step
% slice         struct      Contains all the information about the slice ,e.g Dimension, Position, Index etc
% file          struct      Contains all the loaded file properties ,e.g directory and root name etc
% plin          struct      Contains all the pressure values of the plin slice
% pnl           struct      Contains all the pressure values of the pnl slice
% Bubble        struct      Contains all properties of the Scatterers' Cluster
%
% Outputs:
% pnl           struct      Contains all the pressure values of the pnl slice
% AM, version 200723

function [pnl,plin] = xAM_ConPlot(medium,domain,dslice,file,plin,pnl)
file.plot_waves = 'ind';
harmonics = 2;
if strcmp(file.plot_waves , 'ind')
    i_end = 3;
    name = {'xwave' 'rightwave' 'leftwave'};
    pharm = pnl.xwave_data;
else
    i_end = 1;
    pharm = pnl.leftwave_data(1:size(pnl.xwave_data,1),:,:) + pnl.rightwave_data(1:size(pnl.xwave_data,1),:,:) - pnl.xwave_data(1:size(pnl.xwave_data,1),:,:);
    name = {'Xwave' 'Rightwave' 'Leftwave'};
end
%%
for i = 1:i_end
    subplot_txt1 = 'X Pressure Field';
    if i==2; pharm = pnl.rightwave_data; subplot_txt1 ='Right Wave Pressure Field'; end
    if i==3; pharm = pnl.leftwave_data; subplot_txt1 ='Left Wave Pressure Field'; end
    
    disp(' ')
    disp(['Calculating the ', subplot_txt1,'...'])
    
    subplot_txt =[ subplot_txt1, ' at ',num2str(dslice.savedim(dslice.num)), ' = ', num2str(dslice.pos(dslice.num)) , ' [mm]'];
    dims = size(pharm); dims = dims(2:3);
    plin.fund   = nan(dims);
    pnl.fund   = nan(dims);
    for k=1:harmonics
        eval(['pnl.harm_' name{i}  num2str(k) '= nan(dims);']);
    end
        order = 4;
    %     for iz=1:dims(2)
    %         plin.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
    %         plin.nd((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,1.7*medium.freq0,2.3*medium.freq0,order);
    %         pnl.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert( squeeze(pnl.xwave_data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
    %         for k = 1: harmonics
    %             eval(['pnl.harm_' name{i}  num2str(k) '((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=(maxhilbert(squeeze(pharm(:,:,iz)),domain.dtpar,' num2str(k-0.3) '*medium.freq0,' num2str(k+0.3) '*medium.freq0,order));']);
    %         end
    %     end
    plin.fund = maxhilbert_vec(plin.data,domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
    pnl.fund  = maxhilbert_vec(pnl.xwave_data,domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
    for k = 1: harmonics
        eval(['pnl.harm_' name{i}  num2str(k) '=maxhilbert_vec(pharm,domain.dtpar,' num2str(k-0.3) '*medium.freq0,' num2str(k+0.3) '*medium.freq0,order);']);
        eval(['pnl.diff_' name{i}  num2str(k) '=maxhilbert_vec(20*log10(abs(pnl.xwave_data-plin.data)),domain.dtpar,' num2str(k-0.5) '*medium.freq0,' num2str(k+0.5) '*medium.freq0,order);']);
    end
    %% PLOT THE WAVEFORM SHAPE OF THE 3 DIFFERENT CASES (X,R,L) @ SPECIFIC DEPTH
%     p_interest = pnl.contrastdata;
    disp(['Ploting the signal and frequency spectrum at specific depth ...'])
    colors = { '#0072BD';'#D95319';	'#4DBEEE';	'#77AC30';'#7E2F8E';'#EDB120';'#A2142F';'#009999'};
    
    length_init = 1;
    length_t = 450;
    depth_pos = 148;
    Oversampling_Factor = 5;
    t= linspace(domain.tpar(length_init),domain.tpar(length_t),Oversampling_Factor*(length_t-length_init+1));
    M_I = [1:50:300 301:20:360 361:470 471:20:500 501:50:length(t)];
    f0 = figure('WindowState','maximized');
    hold on; 
    plot(t*1E6,resampleSINC(p_interest        (length_init:length_t,round(dims(1)/2)+1,depth_pos)',Oversampling_Factor)*1E-3,'LineWidth',3)
%     plot(t*1E6,resampleSINC(pnl.rightwave_data(length_init:length_t,round(dims(1)/2)+1,depth_pos)',Oversampling_Factor)*1E-3,'--s','MarkerIndices',M_I,'MarkerSize',4,'LineWidth',3)
%     plot(t*1E6,resampleSINC(pnl.leftwave_data (length_init:length_t,round(dims(1)/2)+1,depth_pos)',Oversampling_Factor)*1E-3,'-.','LineWidth',3)
    set(gca,'FontSize',20)
    xlabel('Time [\mus]')
    ylabel('Acoustic Pressure [kPa]')
    title(['Pressure graph for ' num2str(domain.angle) ' deg at Z = ', num2str(domain.par{3}(depth_pos)), ' [mm]'])
    legend('Both Apertures', 'Right Aperture','Left Aperture')
    grid on;
    %     xlim([0.27 0.55])
    hold off;
    set(gcf,'Color','white')
    set(gca,'FontSize',20)
    %% PLOT THE FREQUENCY SPECTRUM OF THE 3 DIFFERENT CASES (X,R,L) @ SPECIFIC DEPTH
    f02= figure('WindowState','maximized');
    ax = gca;   hold on; grid on;box on;
    p_interest = pnl.contrastdata;
%     pnl.leftwave_data = plin.data;
%     pnl.rightwave_data = 0*plin.data;
    
    F_s = 1/(2*domain.dtpar);
    % Calculate the frequency spectrum based on a sampling frequency and the pressure values
    Nfft_p = 2^nextpow2(2*length_t);                      % Zero-padding to improve fft performance until length equal to 2N
    Y_x = abs(fft(p_interest(length_init:length_t,round(dims(1)/2)+1,depth_pos),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
    f_p = linspace(0,1,Nfft_p/2+1)*F_s;                                                           % [Hz] Frequency sample Number
    % Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
    Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
    S_x = 20*log10(Y_x(Np));
%     
%     Y_left = abs(fft(pnl.leftwave_data(length_init:length_t,round(dims(1)/2)+1,depth_pos),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
%     % Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
%     S_left = 20*log10(Y_left(Np));
%     
%     Y_diff = abs(fft(pnl.xwave_data(length_init:length_t,round(dims(1)/2)+1,depth_pos) - pnl.leftwave_data(length_init:length_t,round(dims(1)/2)+1,depth_pos) - pnl.rightwave_data(length_init:length_t,round(dims(1)/2)+1,depth_pos),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
%     % Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
%     S_diff = 20*log10(Y_diff(Np));
    
    x_plot=plot(f_p*1E-6,S_x - max(S_x),'b-','LineWidth',2);
%     left_plot=plot(f_p*1E-6,S_left - max(S_x),'k--','LineWidth',2);
%     diff_plot=plot(f_p*1E-6,S_diff-max(S_x),'g-.','LineWidth',2);
    legend([x_plot left_plot diff_plot], 'X-wave', 'Right-aperture wave',['Diff @', num2str(domain.par{3}(depth_pos)), ' [mm]' ])
    ylim([-80 0]);
    ylabel(' amplitude [dB]');
    title('Frequency spectrum (INCS)');
    xlabel('Frequency [MHz]');
    xlim([0 F_s*1E-6]);
    set(gca,'FontSize',20);
    set(gcf,'Color','white')
    hold off;
    
    %%  PLOT THE FREQUENCY SPECTRUM OF THE X-WAVE @ DIFFERENT DEPTHS
    disp(['Ploting the frequency spectrum at different depths ...'])
    f02= figure('WindowState','maximized');
    my_col = repmat(linspace(1,0,11)',1,3);  % create your own gray colors color map
    F_s = 1/(2*domain.dtpar);
    % Calculate the frequency spectrum based on a sampling frequency and the pressure values
    Nfft_p = 2^nextpow2(2*length_t);                      % Zero-padding to improve fft performance until length equal to 2N
    f_p = linspace(0,1,Nfft_p/2+1)*F_s;                                                           % [Hz] Frequency sample Number
    Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
    
    depth_of_int = [50:50:domain.dimlen(3)];
    for j = 1:length(depth_of_int)
        length_t = 100;
        if (j<=3) length_t=50;end
        depth_pos = depth_of_int(j);
        Oversampling_Factor = 5;
        t= linspace(domain.tpar(1),domain.tpar(length_t),Oversampling_Factor*length_t);
        ax = gca;   hold on; grid on;box on;
        
        Y_diff = abs(fft(pnl.xwave_data(1:length_t,round(dims(1)/2)+1,depth_pos) - pnl.leftwave_data(1:length_t,round(dims(1)/2)+1,depth_pos) - pnl.rightwave_data(1:length_t,round(dims(1)/2)+1,depth_pos),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
        % Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
        S_diff = 20*log10(Y_diff(Np));
        
        diff_plot(j)=plot(f_p*1E-6,S_diff-max(S_x),'--','LineWidth',2);
        set(diff_plot(j),'Color',my_col(j,:))
        keep_legend{j}=['Diff @', num2str(domain.par{3}(depth_pos)), ' [mm]' ];
    end
    ylim([-80  max(S_diff-max(S_x))]);
    ylabel(' amplitude [dB]');
    title('Frequency spectrum (INCS)');
    xlabel('Frequency [MHz]');
    legend([flip(diff_plot) ] ,flip(keep_legend),'NumColumns',1)
    xlim([0 F_s*1E-6]);
    set(gca,'FontSize',20);
    hold off;
    %% PLOT THE X AMPLITUDE MODULATION [DIFFERENCE IN PRESSURE OF INTEREST IN THE 1H (PLIN,PNL)] @ SPECIFIC DEPTH
    disp(['Ploting the xAM peak pressure through depth...'])
    
    f1 = figure('WindowState','maximized');
    set(gcf, 'Color', 'w');
    hold on ;    grid on;
    
    plot(domain.par{domain.dimval(3)},pnl.fund(round(dims(1)/2)+1,:)*1E-3,'--','Color','black','LineWidth',2)
    plot(domain.par{domain.dimval(3)},plin.fund(round(dims(1)/2)+1,:)*1E-3,'--','Color',[0.8 0.8 0.8],'LineWidth',2)
    %         plot(domain.par{domain.dimval(3)},pnl.harm_xwave1(round(dims(1)/2)+21,:)*1E-3,'--','Color',[0.7 0.7 0.7],'LineWidth',2)
    
    %     for k = 1: 1
    %         eval(['p=plot(domain.par{domain.dimval(3)},pnl.harm_' name{i}  num2str(k) '(round(dims(1)/2)+1,:)*1E-3);']);
    % %         set(p,'Color',sscanf(colors{i}(2:end),'%2x%2x%2x',[1 3])/255,'LineWidth',1.5)
    %     end
    
    xlabel('Depth [mm]')
    ylabel('Acoustic Pressure [kPa] ')
    legend('INCS 20%','INCS 60%','INCS 50%','Location', 'southeast')
    title(['xAM Residual Pressure at ' ,dslice.ylabel(1),' = ',num2str(round(domain.par{domain.dimval(2)}(round(dims(1)/2)+1)*1e+3)/1e+3), ' [mm], ' num2str(domain.angle) ' deg'])
    %     xlim([domain.par{domain.dimval(3)}(1) domain.par{domain.dimval(3)}(end)])
    set(gca,'FontSize',20)
    %%
    disp(['Ploting the xAM percentage ratio through depth...'])
    
    f2 = figure('WindowState','maximized');
    set(gcf, 'Color', 'w');
    hold on ;    grid on;
    
    plot(domain.par{domain.dimval(3)},plin.fund(round(dims(1)/2)+1,:)./pnl.fund(round(dims(1)/2)+1,:)*100,'--','LineWidth',2)
    
    %     xlabel('Depth [mm]')
    %     ylabel('Ratio xAM/x-Wave [%] ')
    %     legend('INCS 70%','INCS 60%','INCS 50%','Location', 'southeast')
    %     title(['xAM Residual Pressure at ' ,dslice.ylabel(1),' = ',num2str(round(domain.par{domain.dimval(2)}(round(dims(1)/2)+1)*1e+3)/1e+3), ' [mm], ' num2str(domain.angle) ' deg'])
    %     xlim([domain.par{domain.dimval(3)}(1) domain.par{domain.dimval(3)}(end)])
    %     set(gca,'FontSize',20)
    %% PLOT THE PRESSURE OF INTEREST IN THE 1H (PLIN,PNL) @ DIFFERENT DEPTHS  @ LATERAL PLANE
    p_interest = squeeze((max(abs(pnl.xwave_data))));
    disp(['Ploting the signal and frequency spectrum at specific depth ...'])
    colors = { '#0072BD';'#D95319';	'#4DBEEE';	'#77AC30';'#7E2F8E';'#EDB120';'#A2142F';'#009999'};
    
    diff_from_mid = 10;
    mid = round(size(p_interest,1)/2+1);
    depth_of_int = [300:20:600];
    Oversampling_Factor = 5;
    x_OS= linspace(domain.par{1}(mid-diff_from_mid),domain.par{1}(mid+diff_from_mid),Oversampling_Factor*(2*diff_from_mid+1));
    f0 = figure('WindowState','maximized');
    hold on;
    clear FWHM_c keep_legend
    max_p_all = max(max(p_interest(mid-diff_from_mid:mid+diff_from_mid,:)));
    for j = 1:length(depth_of_int)
        depth_pos = depth_of_int(j);
        RS_pressure = resampleSINC((p_interest(mid-diff_from_mid:mid+diff_from_mid,depth_pos))',Oversampling_Factor);
        %         RS_pressure = interp1(domain.par{1}(mid-diff_from_mid:mid+diff_from_mid),(p_interest(mid-diff_from_mid:mid+diff_from_mid,depth_pos))',x_OS);
        
        max_p = max(p_interest(mid-diff_from_mid:mid+diff_from_mid,depth_pos));
        [~,maxp_idx] = min(abs(RS_pressure- max_p));
        
        half_max_p = max_p/2;
        [closest_minidx_to_maxp,~]=min(abs(find(RS_pressure-half_max_p<0)- maxp_idx));
        [~,halfmaxp_idx]=min(abs(RS_pressure(max([maxp_idx-closest_minidx_to_maxp-Oversampling_Factor 0]):min([maxp_idx+closest_minidx_to_maxp+Oversampling_Factor length(RS_pressure)])) -half_max_p));
        halfmaxp_idx = halfmaxp_idx+maxp_idx-closest_minidx_to_maxp-Oversampling_Factor-1;
        %         [~,halfmaxp_idx]=min(abs(RS_pressure -half_max_p));
        
        FWHM_c(j,1:2) = [domain.par{3}(depth_pos) ; abs(x_OS(halfmaxp_idx)*2)/medium.dLambdaNN*sind(domain.angle)];
        %         plot(x_OS,(RS_pressure)*1E-3,'LineWidth',2)
        plot(x_OS,20*log10(RS_pressure/max_p),'LineWidth',2)
        keep_legend{j}=['Diff @', num2str(domain.par{3}(depth_pos)), ' [mm]' ];
    end
    
    set(gca,'FontSize',20)
    set(gcf,'Color','white')
    xlabel('X [mm]')
    ylabel('Normalized Amplitude [dB]')
    title(['Pressure, Lateral Profile for ' num2str(domain.angle) ' deg'])
    plot([x_OS(1) x_OS(end)],20*log10(half_max_p*ones(2,1)*1E-3/(max_p*1E-3)),'k--','LineWidth',2)
    legend([keep_legend, 'Half Pressure Amplitude'])
    legend(keep_legend)
    grid on;
    hold off;
    %% ====================================================== Harmonic Plot =================================================
    Nsubplot = 1;
    splot_div = 1;
    Add = 0;
    dBVALUES=3;
    
    FontSize = 20;
    
    %     f2=figure('WindowState','maximized');
    %     imagesc(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},plin.fund)
    %     hold on;
    %     caxis([max(max(plin.fund))-dBVALUES max(max(plin.fund))]+Add)
    %     title('Linear Pressure Field')
    %     xlabel(dslice.xlabel )
    %     ylabel(dslice.ylabel )
    %     BubbleCluster_Colormaps(file)
    %     c = colorbar;
    %     c.Label.String = 'Pressure [MPa]'; c.Label.FontSize = FontSize;
    %     set(c,'Ticklabels',get(c,'Ticklabels'))
    %     set(gca,'YDir','normal')
    %     if dslice.savedim(dslice.num) =='z'; set(gca,'XDir','reverse');set(gca,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
    %     set(gcf, 'Color', 'white');
    %     set(gca,'FontSize',FontSize);
    %     hold off;
    
    for k = 1: harmonics
        if (Nsubplot~=1 && k-4*floor((k-1)/4)==1); f3=figure('WindowState','maximized'); subplot(Nsubplot/splot_div,Nsubplot/(Nsubplot/splot_div),k-4*floor((k-1)/4)); ...
                sgtitle(subplot_txt,'FontSize',FontSize+3);  spl = lower(split(subplot_txt1,' ')); ...
        else; figure('WindowState','maximized');end
    
%     pnl_harmi = eval(['pnl.harm_' name{i}  num2str(k)])*1E-3;
    pnl_harmi= eval(['pnl.diff_' name{i}  num2str(k)]);
%     pnl_harmi = 20*log10(squeeze(max(abs(plin.data))));
    imagesc(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},pnl_harmi)
    
    hold on;
    caxis([max(max(pnl_harmi))-dBVALUES max(max(pnl_harmi))]+Add)
%     title(['Beam Profile - Pressure, ' num2str(domain.angle) ' deg, ' name{i} ', ' num2str(k) 'H, 60% Butterworth Filter'])
    xlabel(dslice.xlabel )
    ylabel(dslice.ylabel )
    xAM_Colormaps(file)
    c = colorbar;
    c.Label.String = 'Amplitude [dB]'; c.Label.FontSize = FontSize;
    set(gca,'YDir','normal')
    if dslice.savedim(dslice.num) =='z'; set(gca,'XDir','reverse');set(gca,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
    
    ax = gca;
    ax.FontSize = FontSize;
    set(gcf, 'Color', 'white');
    hold off;
    
    end
    
    if (strcmp(file.saveplot,'yes'))
        figure(f1)
        folder1 = ['../',file.savedir,'/',file.dirname,'_plot_',spl{1},'_',spl{2},'_slice',int2str(dslice.num),'.png'];
        cd export_fig
        export_fig  img1 -painters -q110
        cd ../
        figure(f2)
        folder2 = ['../',file.savedir,'/',file.dirname,'_linear_',spl{2},'_slice',int2str(dslice.num),'.png'];
        cd export_fig
        export_fig  img2 -painters -q110
        cd ../
        figure(f3)
        folder3 = ['../',file.savedir,'/',file.dirname,'_',spl{1},'_',spl{2},'_slice',int2str(dslice.num),'.png'];
        cd export_fig/
        export_fig  img3 -painters -q110
        movefile('img1.png', folder1)
        movefile('img2.png', folder2)
        movefile('img3.png', folder3)
        cd ../
    end
    
    
end
end