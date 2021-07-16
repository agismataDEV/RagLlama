% CREATING MATRICES FOR NONLINEAR AND CONTRAST DATA AND PLOTTING THE PRESSURE FIELD
%
% [pnl] = BubbleCluster_ConPlot(medium,domain,slice,file,plin,pnl,Bubble)
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

function [pnl,plin] = BubbleCluster_ConPlot(medium,domain,dslice,file,plin,pnl,Bubble)

if strcmp(file.plot_contrast,'yes')
    i_end = 2;
else
    i_end = 1;
end


pharm = pnl.data;
%%
Nsubplot = 4;
splot_div =2;
dBVALUES=20;
Add = 0;
harmonics = domain.Fnyq -1 ;
for i = 1:i_end
    subplot_txt1 = 'Total Pressure Field';
    if i==2; pharm = pnl.contrastdata; subplot_txt1 ='Contrast Pressure Field'; end
    
    disp(' ')
    disp(['Calculating the ', subplot_txt1,'...'])
    
    subplot_txt =[ subplot_txt1, ' at ',num2str(dslice.savedim(dslice.num)), ' = ', num2str(dslice.pos(dslice.num)) , ' [mm]'];
    dims = size(plin.data); dims = dims(2:3);
    plin.fund   = nan(dims);
    pnl.atten   = nan(dims);
    plin.atten  = nan(dims);
    for i=1:harmonics
        eval(['pnl.harm' num2str(i) '= nan(dims);']);
    end
    
    order = 4;
    % There is not significant change if I use the 1st block of code for a cluster of bubbles,
    % But the opposite cannot happen ,because the 2nd block of code can not accurately visualize
    % the behaviour of one bubble , due to the spectral leakage in the lower harmonics
    % generate fundamental, higher harmonic and superharmonic profiles
    if Bubble.N>0
        
        % for iz=1:dims(2)
            % plin.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=20*log10(maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order));
            % for i = 1: harmonics
                % eval(['pnl.harm' num2str(i) '((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=20*log10(maxhilbert(squeeze(pharm(:,:,iz)),domain.dtpar,' num2str(i-0.3) '*medium.freq0,' num2str(i+0.3) '*medium.freq0,order));']);
            % end
                        % pnl.atten((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(pnl.data(:,:,iz)),domain.dtpar,.1*medium.freq0,(domain.Fnyq-0.3)*medium.freq0,order);
                        % plin.atten((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.1*medium.freq0,(domain.Fnyq-0.3)*medium.freq0,order);
                        % pnl.atten_fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(pnl.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
                        % plin.atten_fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
            
        % end
		plin.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),:)=maxhilbert_vec(plin.data,domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
		plin.nd((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),:)  =maxhilbert_vec(plin.data,domain.dtpar,1.7*medium.freq0,2.3*medium.freq0,order);
		pnl.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),:) =maxhilbert_vec( pnl.xwave_data,domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
		for i = 1: harmonics
            eval(['pnl.harm' num2str(i) '((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),:)=20*log10(maxhilbert_vec(pharm,domain.dtpar,' num2str(i-0.3) '*medium.freq0,' num2str(i+0.3) '*medium.freq0,order));']);
        end
        
    end
    
    %% PLOTS
    disp(['Ploting the ', subplot_txt1,'...'])
    
    FontSize = 20;
    %     I = 20*log10(pnl.atten./plin.atten);
    %     pnl.attenuation(dslice.num,:) = [dslice.pos(dslice.num); mean(mean(I))/(dslice.pos(dslice.num) - min(Bubble.LocGlob(:,3)))*10];  % [mm] to [cm]
    %     I_fund = 20*log10(pnl.atten_fund./plin.atten_fund);
    %     pnl.attenuation_fund(dslice.num,:) = [dslice.pos(dslice.num); mean(mean(I_fund))/(dslice.pos(dslice.num) - min(Bubble.LocGlob(:,3)))*10];  % [mm] to [cm]
    colors = { '#0072BD';'#D95319';	'#4DBEEE';	'#77AC30';'#7E2F8E';'#EDB120';'#A2142F';'#009999'};

    f1 = figure('WindowState','maximized');
    set(gcf, 'Color', 'w');
    hold on
    plot(domain.par{domain.dimval(3)},plin.fund(round(dims(1)/2)+1,:),'--','Color','black')
    for i = 1: harmonics
        eval(['p=plot(domain.par{domain.dimval(3)},pnl.harm' num2str(i) '(round(dims(1)/2)+1,:));']);
        set(p,'Color',sscanf(colors{i}(2:end),'%2x%2x%2x',[1 3])/255,'LineWidth',1.5)
    end
    
    grid on
    xlabel(dslice.xlabel)
    ylabel('Pressure [dB re 1 Pa] ')
    legendtxt = {'Linear';'NonLin Fund'};
    for i =2:harmonics 
        legendtxt{i+1} = eval(['"',num2str(i),'H";']);
    end
    
    legend(legendtxt,'Location', 'southeast','Orientation','horizontal','NumColumns',3)
    titletxt=['Pressure at ' ,dslice.ylabel(1),' = ',num2str(round(domain.par{domain.dimval(2)}(round(dims(1)/2)+1)*1e+3)/1e+3), ' [mm], at slice ', upper(num2str(dslice.savedim(dslice.num))),' = ',num2str(round(dslice.pos(dslice.num)*10)/10),' [mm]'];
    title(titletxt)
    xlim([domain.par{domain.dimval(3)}(1) domain.par{domain.dimval(3)}(end)])
    ylim([0 max(max(plin.fund))]);
        set(gcf, 'Color', 'white');
    set(gca,'FontSize',FontSize);
    
    %====================================================== Harmonic Plot =================================================
    Nsubplot = 1;
    splot_div =1;
    dBVALUES=30;
    Add = 0;
    
    f2=figure('WindowState','maximized');
    imagesc(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},plin.fund)
    hold on;
    caxis([max(max(plin.fund))-dBVALUES max(max(plin.fund))]+Add)
    title(['Linear Pressure Field'])
    xlabel(dslice.xlabel )
    ylabel(dslice.ylabel )
    BubbleCluster_Colormaps(file)
    c = colorbar;
    c.Label.String = 'Pressure [dB]'; c.Label.FontSize = FontSize;
    set(gca,'YDir','normal')
    if dslice.savedim(dslice.num) =='z'; set(gca,'XDir','reverse');set(gca,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
    set(gcf, 'Color', 'white');
    set(gca,'FontSize',FontSize);
    hold off;
    
    for i = 1: 1
        if (i-4*floor((i-1)/4)==1); f3=figure('WindowState','maximized');end
        pnl_harmi = eval(['pnl.harm' num2str(i)]);
        pnl_harmi = pnl.FS;
        subplot(Nsubplot/splot_div,Nsubplot/(Nsubplot/splot_div),i-4*floor((i-1)/4))
        imagesc(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},pnl_harmi)
        
        hold on;
        caxis([max(max(pnl_harmi))-dBVALUES max(max(pnl_harmi))]+Add)
        if Bubble.N >6
            rectangle('Position',[domain.min_dim(domain.dimval(3)) domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(3))-domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(2))-domain.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor','white','LineWidth',2)
        elseif Bubble.N>=1
            plot(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),'xk','MarkerSize',10,'LineWidth',3)
        end
        title([num2str(i) 'H'])
        xlabel(dslice.xlabel )
        ylabel(dslice.ylabel )
        BubbleCluster_Colormaps(file)
        c = colorbar;
        c.Label.String = 'Pressure [dB]'; c.Label.FontSize = FontSize;
        set(gca,'YDir','normal')
        if dslice.savedim(dslice.num) =='z'; set(gca,'XDir','reverse');set(gca,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
        
        ax = gca;
        ax.FontSize = FontSize;
        set(gcf, 'Color', 'white');
        sgtitle(subplot_txt,'FontSize',FontSize+3);
        spl = lower(split(subplot_txt1,' '));
        hold off;
        
    end
    
    if (strcmp(file.saveplot,'yes'))
        figure(f1)
        folder1 = ['../',file.savedir,'/',file.dirname,'_plot_',spl{1},'_',spl{2},'_slice',int2str(dslice.num),'.png'];
        cd export_fig
        export_fig  img1 -q110
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