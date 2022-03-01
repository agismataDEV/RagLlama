% ***********************************************************************************************************
% PLOTING THE LOCATION OF EACH SCATTERER INDIVIDUALLY AND OF THE CLUSTER AS A WHOLE
%
% **********************************************************
% [domain] = BubbleCluster_LocPlot(domain,Bubble,slice,file)
% **********************************************************
% Plot the domain with the existence of the scatterers
% Plot all the sideviews of the domain to understand exactly where the scatterers
% are placed
%
% Inputs:
% domain        struct      Contains all the computational domain properties, e.g space and time step
% slice         struct      Contains all the information about the slice ,e.g Dimension, Position, Index etc
% file          struct      Contains all the loaded file properties ,e.g directory and root name etc
% Bubble        struct      Contains all properties of the Scatterers' Cluster
%
% Outputs:
% domain        struct      Contains all the computational domain properties, e.g space and time step
%
% ******************
% AM, version 200723
% ******************
%
% ***********************************************************************************************************

function [domain] = BubbleCluster_LocPlot(domain,dslice,file,Bubble,name,MarkerStyle)

%%
disp('Ploting the microbubble cluster...')
MarkerSize =  8/(log10(Bubble.N)+1);
if (Bubble.N >0 )
    %=============================== 3D Domain View ==============================
    f4 = figure('WindowState','maximized');
    LocPlot_FontSize = 15;
    subplot(2,3,1)
    hold on;
    box on;
    [Z,X] = meshgrid(domain.par{3},domain.par{1});
    s=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],dslice.pos(dslice.num)*ones(4,4) ,...
        'EdgeColor','black', 'FaceAlpha' ,0.6);
    % shading interp;
    if (domain.dimval(1)==1)
        rotate(s,[1 0 0 ], 90)
    elseif (domain.dimval(1)==2)
        rotate(s,[0 0 1 ], 0)
    elseif (domain.dimval(1)==3)
        rotate(s,[0 1 0], -90)
    end
    %     [ y x z]
    p=plot3(Bubble.LocGlob(:,3),Bubble.LocGlob(:,1),Bubble.LocGlob(:,2),MarkerStyle,'MarkerSize',MarkerSize);
    if (strcmp(file.plot_colour,'gray')); s.FaceColor =[ 192 192 192]/255; p.MarkerEdgeColor = sscanf('404040','%2x%2x%2x',[1 3])/255 ;  end
    
    % Simplified visualization of transducer with 8 elements
    El = 8;
    y_trans(1,1)=domain.par{1}(1);
    y_trans(1,2)=domain.par{1}(1) + (domain.par{1}(end)-domain.par{1}(1))/(El*1.2);
    for i = 1:El
        plot3(Z(1,1)*ones(5,1) ,[y_trans(i,1) y_trans(i,1) y_trans(i,2) y_trans(i,2) y_trans(i,1)]', ...
            [domain.par{2}(1) domain.par{2}(end) domain.par{2}(end) domain.par{2}(1) domain.par{2}(1)]'+ 0.1*(domain.par{2}(end)-domain.par{2}(1))*[1 -1 -1 1 1]', 'k')
        y_trans(i+1,1) = y_trans(i,2) +0.2*(domain.par{1}(end)-domain.par{1}(1))/(El*1.2);
        y_trans(i+1,2) = y_trans(i+1,1) + (domain.par{1}(end)-domain.par{1}(1))/(El*1.2);
    end
    view(-37.5,30)
    
    xlim([min(domain.par{3}) max(domain.par{3})])
    ylim([min(domain.par{1}) max(domain.par{1})])
    zlim([min(domain.par{2}) max(domain.par{2})])
    
    x = xlabel('Z [mm]');
    set(x, 'Position',get(x, 'Position').*[1,1,1],'Rotation',10)
    y = ylabel('X [mm]');
    set(y, 'Position',get(y, 'Position').*[1,1,1],'Rotation',-21)
    zlabel('Y [mm]')
    title('Total Domain View')
    grid on
    hold off;
    
    set(gca,'FontSize',LocPlot_FontSize)
    %======================= 3D Microbubble Cluster Focused position =============
    subplot(2,3,2)
    hold on;
    box on;
    [Z,X] = meshgrid(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)});
    s=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],dslice.pos(dslice.num)*ones(4,4),...
        'EdgeColor','black', 'FaceAlpha' ,0.6);
    p = plot3(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),Bubble.LocGlob(:,domain.dimval(1)),MarkerStyle,'MarkerSize',MarkerSize );
    if (strcmp(file.plot_colour,'gray')); s.FaceColor =[ 192 192 192]/255; p.MarkerEdgeColor =  sscanf('404040','%2x%2x%2x',[1 3])/255 ;end
    view(-37.5,30)
    
    xlim(round(10*([min(Bubble.LocGlob(:,domain.dimval(3))) max(Bubble.LocGlob(:,domain.dimval(3)))] + [-domain.dzpar domain.dzpar]))/10)
    ylim(round(10*([min(Bubble.LocGlob(:,domain.dimval(2))) max(Bubble.LocGlob(:,domain.dimval(2)))] + [-domain.dxpar domain.dxpar]))/10)
    zlim(round(10*([min(Bubble.LocGlob(:,domain.dimval(1))) max(Bubble.LocGlob(:,domain.dimval(1)))] + [-domain.dypar domain.dypar]))/10)
    % xlim([ 40 45])
    % ylim([ -4 -3])
    % zlim([ 1 2])
    
    x = xlabel(dslice.xlabel);
    set(x, 'Position',get(x, 'Position').*[1,1,1],'Rotation',10)
    y = ylabel(dslice.ylabel);
    set(y, 'Position',get(y, 'Position').*[1,1,1],'Rotation',-21)
    zlabel(dslice.zlabel)
    title([name,' Cluster Focused Region'])
    grid on
    hold off;
    
    set(gca,'FontSize',LocPlot_FontSize)
    %====================== X-Z Side View of Focused Region ==========================
    subplot(2,3,3)
    
    hold on;
    p=plot(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(1)),MarkerStyle,'MarkerSize',MarkerSize);
    
    % Create dimensions that encloses the microbubble cluster but includes the
    % slice of the domain stored in the output file
    domain.min_dim = min([Bubble.LocGlob ; domain.par{1}(end) domain.par{2}(end) domain.par{3}(end)]);
    domain.max_dim = max([Bubble.LocGlob ; domain.par{1}(1) domain.par{2}(1) domain.par{3}(1)]);
    domain.min_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)) = min(domain.min_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)), dslice.pos(dslice.num));
    domain.max_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)) = max(domain.max_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)), dslice.pos(dslice.num));
    
    rect_w =  0.01*(domain.max_dim(domain.dimval(1))+domain.dypar - (domain.min_dim(domain.dimval(1))-domain.dypar));
    r = rectangle('Position',[floor(10*(domain.min_dim(domain.dimval(3))-domain.dzpar))/10 dslice.pos(dslice.num)-rect_w  ceil(10*(domain.max_dim(domain.dimval(3))+domain.dzpar-(domain.min_dim(domain.dimval(3))-domain.dzpar)))/10 rect_w],...
        'EdgeColor','k', 'FaceColor' , [55 90 83]/255);%'#8BE5D3'); %% Visualization of slice
    if (strcmp(file.plot_colour,'gray')); p.MarkerEdgeColor =  sscanf('404040','%2x%2x%2x',[1 3])/255 ; r.FaceColor = [ 192 192 192]/255;end
    
    xlim(round(10*([ domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(3))] + [-domain.dzpar domain.dzpar]))/10)
    ylim(round(10*([ domain.min_dim(domain.dimval(1)) domain.max_dim(domain.dimval(1))] + [-domain.dypar domain.dypar]))/10)
    
    xlabel(dslice.xlabel)
    ylabel(dslice.zlabel)
    
    title([dslice.zlabel(1),' - ',dslice.xlabel(1), ' View of Focused Region'])
    grid on
    hold off;
    
    set(gca,'FontSize',LocPlot_FontSize)
    %====================== 2D Side View of Focused Region ==========================
    subplot(2,3,[4 5])
    hold on;
    p = plot(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),MarkerStyle,'MarkerSize',MarkerSize );
    if (strcmp(file.plot_colour,'gray')); p.MarkerEdgeColor =  sscanf('404040','%2x%2x%2x',[1 3])/255 ;end
    
    xlim(round(10*[min(domain.par{domain.dimval(3)}) max(domain.par{domain.dimval(3)})])/10)
    ylim(round(10*[min(domain.par{domain.dimval(2)}) max(domain.par{domain.dimval(2)})])/10)
    
    xlabel(dslice.xlabel)
    ylabel(dslice.ylabel);
    
    title([dslice.ylabel(1),' - ',dslice.xlabel(1), ' View of Focused Region'])
    grid on
    hold off;
    
    set(gca,'FontSize',LocPlot_FontSize)
    %====================== Z-Y Side View of Focused Region ==========================
    subplot(2,3,6)
    hold on;
    p = plot(Bubble.LocGlob(:,domain.dimval(2)),Bubble.LocGlob(:,domain.dimval(1)),MarkerStyle,'MarkerSize',MarkerSize );
    
    r = rectangle('Position',[floor(10*(domain.min_dim(domain.dimval(2))-domain.dxpar))/10 dslice.pos(dslice.num)-rect_w  ...
        ceil(10*(domain.max_dim(domain.dimval(2))+domain.dxpar-(domain.min_dim(domain.dimval(2))-domain.dxpar)+domain.dxpar))/10 rect_w],...
        'EdgeColor','k', 'FaceColor'  , [55 90 83]/255);%'#8BE5D3') ;%% Visualization of slice
    if (strcmp(file.plot_colour,'gray')); p.MarkerEdgeColor =  sscanf('404040','%2x%2x%2x',[1 3])/255 ; r.FaceColor = [ 192 192 192]/255;end
    
    xlim(round(10*([ domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(2))] + [-domain.dxpar domain.dxpar]))/10)
    ylim(round(10*([ domain.min_dim(domain.dimval(1)) domain.max_dim(domain.dimval(1))] + [-domain.dypar domain.dypar]))/10)
    % xlim([-0.1 0.1])
    % ylim([-0.1 0.1])
    
    xlabel(dslice.ylabel);
    ylabel(dslice.zlabel);
    
    title([dslice.zlabel(1),' - ',dslice.ylabel(1), ' View of Focused Region'])
    grid on
    hold off;
    
    %============================= Create subplot title =======================
    integerPart = num2str(Bubble.N);
    % Reverse the integer-part string.
    integerPart=integerPart(end:-1:1);
    % Insert commas every third entry.
    integerPart=[sscanf(integerPart,'%c',[3,inf])' repmat(',',ceil(length(integerPart)/3),1)]';
    integerPart=integerPart(:)';
    % Strip off any trailing commas.
    integerPart=deblank(integerPart(1:(end-1)));
    integerPart = integerPart(end:-1:1);
    
    sgtitle([name,' Cluster of ',integerPart,' Bubble(s) per ',num2str(round(Bubble.Volume*1e-3,2)),' [mL] at ',num2str(dslice.savedim(dslice.num)), ' = ', num2str(dslice.pos(dslice.num)) , ' [mm]'],'FontSize',LocPlot_FontSize+5);
    
    % set(gcf, 'Color', 'w');
    set(gca, 'Color', 'none'); % Sets axes background
    set(gca,'FontSize',LocPlot_FontSize)
    %========================== Save plot ===================================
    % if (strcmp(file.saveplot,'yes'))
    %     saveas(f4,[file.savedir,'/',file.dirname,'_cluster_slice',int2str(dslice.num),'.jpg'])
    % end
    
    %% Microbubble Radius distribution
    %%
    if (sum(Bubble.Radius==0)==0)
        f5 = figure('WindowState','maximized');
        BinNum  = 100;
        histfit(Bubble.Radius,BinNum,'Gamma');
        h = get(gca,'Children');
        set(h(2),'FaceColor',[.8 .8 .8]) % Bar Properties
        set(h(1),'Color','blue')         % Line Properties
        %     set(h,'FaceColor','#404040')
        ax = gca;
        xlim([min(Bubble.Radius)-1E-10 max(Bubble.Radius)])
        set(gca,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Radius [\mum]';
        ax.YLabel.String = 'No. of Microbubbles';
        title('Microbubble Radius Distribution')
        set(gca, 'Color', 'white'); % Sets axes background
        set(gcf, 'Color', 'white'); % Sets axes background
        set(gca,'FontSize',30)
    end
    
    %%
    if (strcmp(file.saveplot,'yes'))
        figure(f4)
        folder = ['../',file.savedir,'/',file.dirname,'_cluster_slice',int2str(dslice.num),'.jpg'];
        cd export_fig/
        export_fig  folder -painters -q110 -transparent
        movefile('folder.png',folder)
        cd ../
    end
end
end

