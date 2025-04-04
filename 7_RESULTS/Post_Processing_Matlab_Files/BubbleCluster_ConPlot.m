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
harmonics = 4 ;
dims = size(plin.data); dims = dims(2:3);
%%
for i = 1:i_end
    subplot_txt1 = 'Total Pressure Field';
    if i==2; pharm = pnl.contrastdata; subplot_txt1 ='Contrast Pressure Field'; end
    
    disp(' ')
    disp(['Calculating the ', subplot_txt1,'...'])
    
    subplot_txt =[ subplot_txt1, ' at ',num2str(dslice.savedim(dslice.num)), ' = ', num2str(dslice.pos(dslice.num)) , ' [mm]'];
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

%         
        for iz=1:dims(2)
            plin.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=20*log10(maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order));
            pnl.fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=20*log10(maxhilbert(squeeze(pnl.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order));
            for i = 1: harmonics
                pnl.harm{i}((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=20*log10(maxhilbert(squeeze(pharm(:,:,iz)),domain.dtpar,(i-0.3)*medium.freq0,(i+0.3)*medium.freq0,order));
            end
%                         pnl.atten((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(pnl.data(:,:,iz)),domain.dtpar,.1*medium.freq0,(domain.Fnyq-0.3)*medium.freq0,order);
%                         plin.atten((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.1*medium.freq0,(domain.Fnyq-0.3)*medium.freq0,order);
%                         pnl.atten_fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(pnl.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
%                         plin.atten_fund((1:dims(1))+domain.TXYoff(iz,2)-min(domain.TXYoff(:,2)),iz)=maxhilbert(squeeze(plin.data(:,:,iz)),domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order);
            
        end
% 		plin.fund =20*log10(maxhilbert_vec(plin.data,domain.dtpar,.7*medium.freq0,1.3*medium.freq0,order));
% 		for i = 1: harmonics
%             pnl.harm{i}=20*log10(maxhilbert_vec(pharm,domain.dtpar,(i-0.3)*medium.freq0,(i+0.3)*medium.freq0,order));
%         end
    
    %% PLOTS
   disp(['Ploting the ', subplot_txt1,'...'])
    
%     I = 20*log10(pnl.atten./plin.atten);
%     pnl.attenuation(dslice.num,:) = [dslice.pos(dslice.num); mean(mean(I))/(dslice.pos(dslice.num) - min(Bubble.LocGlob(:,3)))*10];  % [mm] to [cm]
%     I_fund = 20*log10(pnl.atten_fund./plin.atten_fund);
%     pnl.attenuation_fund(dslice.num,:) = [dslice.pos(dslice.num); mean(mean(I_fund))/(dslice.pos(dslice.num) - min(Bubble.LocGlob(:,3)))*10];  % [mm] to [cm]
    colors = { '#0072BD';'#D95319';	'#4DBEEE';	'#77AC30';'#7E2F8E';'#EDB120';'#A2142F';'#009999'};
    
    f1 = figure('WindowState','maximized');
    set(gcf, 'Color', 'w');
    hold on
    
    plin_plot  = plin.fund;20*log10(squeeze(max(abs(pnl.data))));
    plot(domain.par{domain.dimval(3)},plin_plot(round(dims(1)/2)+1,:),'--','Color','black')
    for i = 1: harmonics
        eval(['p=plot(domain.par{domain.dimval(3)},pnl.harm' num2str(i) '(round(dims(1)/2)+1,:));']);
        set(p,'Color',sscanf(colors{i}(2:end),'%2x%2x%2x',[1 3])/255,'LineWidth',1.5)
    end
%     xline(Bubble.LocRange(3,1),'--','LineWidth',3); xline(Bubble.LocRange(3,2),'--','LineWidth',3)
    grid on
    xlabel(dslice.xlabel)
    ylabel('Pressure [dB re 1 Pa] ')
    legend('Linear P0','NonLin P0', '2H', '3H', '4H','Location', 'southeast')
    title(['Pressure at ' ,dslice.ylabel(1),' = ',num2str(round(domain.par{domain.dimval(2)}(round(dims(1)/2)+1)*1e+3)/1e+3), ' [mm]'])
    xlim([domain.par{domain.dimval(3)}(1) domain.par{domain.dimval(3)}(end)])
 
    %%
    %====================================================== Harmonic Plot =================================================
    Nsubplot = 1;
    splot_div =1;
    dBVALUES=20;
    harmonics = 1;
    Add = 0;
    i=10;
    j=4;
	Bubble=ScPop{1};
    PS=ScPop{2};
%     plin_plot = squeeze(20*log10(max(abs(pnl.data_40{i,j}-plin.data_40{i}))));
    plin_plot = squeeze(20*log10(max(abs(hilbert(plin.data)))));
        FontSize = 25;
    
    f4=figure('WindowState','maximized');
    imagesc(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},plin_plot)
    hold on;
    caxis([max(max(plin_plot))-dBVALUES max(max(plin_plot ))]+Add)
%     if Bubble.N >6
%         rectangle('Position',[Bubble.min_dim(domain.dimval(3)) Bubble.min_dim(domain.dimval(2)) Bubble.max_dim(domain.dimval(3))-Bubble.min_dim(domain.dimval(3)) Bubble.max_dim(domain.dimval(2))-Bubble.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor','white','LineWidth',4)
%     elseif Bubble.N>=1
%         plot(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),'x','Color','white','MarkerSize',10,'LineWidth',3)
%     end
%    if PS.N >6
%         rectangle('Position',[PS.min_dim(domain.dimval(3)) PS.min_dim(domain.dimval(2)) PS.max_dim(domain.dimval(3))-PS.min_dim(domain.dimval(3)) PS.max_dim(domain.dimval(2))-PS.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor','red','LineWidth',2)
%     elseif PS.N>=1
%         plot(PS.LocGlob(:,domain.dimval(3)),PS.LocGlob(:,domain.dimval(2)),'xk','MarkerSize',10,'LineWidth',3)
%     end
    % First MB Cloud
%     r = rectangle('Position',[ScPop{2}.LocRange(3,1) ScPop{2}.LocRange(1,1) diff(ScPop{2}.LocRange(3,:))  diff(ScPop{2}.LocRange(1,:))],...
%         'EdgeColor',[0.64,0.08,0.18], 'FaceColor' , 'none','LineWidth',4,'LineStyle','--');%'#8BE5D3'); %% Visualization of slice
    
    % Second MB Cloud
%     [xunit2,yunit2] = circle([sum(ScPop{1}.LocRange(3,:))/2, sum(ScPop{1}.LocRange(1,:))/2], diff(ScPop{1}.LocRange(3,:))/2);
%     plot(xunit2,yunit2, '--','Color',"#0072BD",'LineWidth',3);
    
    title('Linear Pressure Field')
        title(['F0, [0.7 - 1.3] MHz'] )
        title(['2H, [1.7 - 2.3] MHz'] )
        title(['3H, [2.7 - 3.3] MHz'] )
%         title(['Incident Pressure Field, 1.7MHz, SRC ',num2str(i/2),'$\lambda$'] )
%         title(['Scattered Pressure Field, 1.7MHz, SRC ',num2str(i/2),'$\lambda$, MB ',num2str(j),'$\lambda$'] )
        title({'Scattered Pressure Field,', 'R=4 $\mu$m, f$_{r}$=0.7MHz'} )
    xlabel(dslice.xlabel)
    ylabel(dslice.ylabel )
    BubbleCluster_Colormaps(file)
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String = 'Pressure [dB]'; c.Label.FontSize = FontSize;
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.03))
    set(gca,'YDir','normal')
    if dslice.savedim(dslice.num) =='z'; set(gca,'XDir','reverse');set(gca,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
    set(gcf, 'Color', 'white');
    set(gca,'FontSize',FontSize);
    hold off;
    imgrotate(gca,gcf);    set(gca,'YDir','reverse')
%     caxis([-19 1])
    axis image
    %%
%     caxislim = [57  87 ; 45.0489   70.0489; 40.1915   65.1915; 36.0642   61.0642];
 harmonics=3;
 subplot_txt={'F0 [0.7 ,1.3] MHz','2H [1.7 ,2.3] MHz','3H [2.7 ,3.3] MHz','4H [3.7 ,4.3] MHz'};
 f3=figure('WindowState','maximized');
 t=tiledlayout(ceil(harmonics/4),harmonics/(ceil(harmonics/4)),'TileSpacing','compact','Padding','tight');
    for i = 1: harmonics
        nexttile
        pnl_harmi = pPD.harm{i};   
        imagesc(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},pnl_harmi)
        caxis([max(max(p1.harm{i}))-dBVALUES max(max(p1.harm{i}))]+Add)
        
        hold on;
        if Bubble.N >6
            rectangle('Position',[Bubble.min_dim(domain.dimval(3)) Bubble.min_dim(domain.dimval(2)) Bubble.max_dim(domain.dimval(3))-Bubble.min_dim(domain.dimval(3)) Bubble.max_dim(domain.dimval(2))-Bubble.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor','white','LineWidth',2)
        elseif Bubble.N>=1
            plot(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),'xk','MarkerSize',10,'LineWidth',3)
        end
        title([num2str(i) 'H'] )
        if (i==1) ;title('F0' ); end
        xlabel(dslice.xlabel)
        ylabel(dslice.ylabel)
        BubbleCluster_Colormaps(file)
        c = colorbar;
        c.Label.Interpreter = 'latex';
        c.Label.String = 'Pressure [dB]'; c.Label.FontSize = FontSize;
        set(gca,'YDir','normal')
        if dslice.savedim(dslice.num) =='z'; set(gca,'XDir','reverse');set(gca,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
        
        ax = gca;
        ax.FontSize = FontSize;
        set(gcf, 'Color', 'white');
        title(subplot_txt{i},'FontSize',FontSize+3);
        hold off;
        imgrotate(gca,gcf);    set(gca,'YDir','reverse')
        if (i>1); q=gca;q.YTickLabel=''; q.YLabel.Visible='off';end
%         axis image;
        
    end
     title(t,{'Scattered Pressure Field, R=2.8 $\mu$m, f$_{r}$=1MHz'},'FontSize',25,'Interpreter','latex')
     title(t,{'Scattered Pressure Field, R=[0, 20] $\mu$m'},'FontSize',25,'Interpreter','latex')
    %%
    
    %%
    
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
        cd ../
        figure(f4)
        folder4 = ['../',file.savedir,'/',file.dirname,'_linear_cluster_',spl{2},'_slice',int2str(dslice.num),'.png'];
        cd export_fig/
        export_fig  img4 -painters -q110
        movefile('img1.png', folder1)
        movefile('img2.png', folder2)
        movefile('img3.png', folder3)
        movefile('img4.png', folder4)
        cd ../
    end
    
    
end
end