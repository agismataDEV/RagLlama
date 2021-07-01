% CREATING MATRICES FOR NONLINEAR AND CONTRAST DATA AND PLOTTING THE PRESSURE FIELD
%
% Bubble_PlayMovies(medium,domain,dslice,file,plin,pnl,Bubble)
% Play movies from the extracted results
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
% AM, version 201005

function Bubble_PlayMovies(medium,domain,dslice,file,plin,pnl,Bubble)

if (strcmp(file.play_movies,'yes'))
    mov = VideoWriter('PressureFieldEvolution.avi','Motion JPEG AVI');
    mov.FrameRate = 30; mov.Quality = 75;
    open(mov);
    
    disp(' ')
    disp(['Play movie ...'])
    
    sizeP = size(plin.data);
    Linear_Field=-500*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    Scatter_Field=-500*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    Total_Field=-500*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    
    if (domain.dimval(1) ==1 || domain.dimval(1)==2)   % For X and Y , the domain should be moved in order to get a nice movie
%         Scatter_Field=round(ifft(fft(reshape(pnl.contrastdata,sizeP(1)*sizeP(2),sizeP(3))',[],2).* exp(-2i*pi/10*[1:sizeP(3)]),[],2))
        for k=1:sizeP(3)
                       Linear_Field(k:sizeP(1)+k-1,:,k)=20*log10(abs(plin.data(:,:,k)));
                       Scatter_Field(k:sizeP(1)+k-1,:,k)=20*log10(abs(pnl.contrastdata(:,:,k)));
            Total_Field(k:sizeP(1)+k-1,:,k)=20*log10(abs(pnl.data(:,:,k)));
        end
        
        i_end = sizeP(1)+sizeP(3); % For X And Y the plot is a moving window at each Z position
        
    else   % For Z , the movie can be ploted straigh away
        for i_t=1:sizeP(1)
            Linear_Field(i_t+sizeP(3),:,:) =20*log10(abs(plin.data(i_t,:,:)));%./(max(max(max(plin))))));
            Scatter_Field(i_t+sizeP(3),:,:)=20*log10(abs(pnl.contrastdata(i_t,:,:)));%./(max(max(max(plin))))));
            Total_Field(i_t+sizeP(3),:,:) =20*log10(abs(pnl.data(i_t,:,:)));%./(max(max(max(plin.data))))));
        end
        i_end = sizeP(1)+sizeP(3);   % For Z the plot is only for the time points
    end
    
    dBVALUES = 70;
    plot_value = Scatter_Field;
    
    %% ========================= Initialize first frame=========================
    % This is done to speedup the process of making images , the most time is lost in making the plots
    num_frames = i_end;
    F(num_frames) = struct('cdata',[],'colormap',[]);
    
    f2 = figure('WindowState','maximized');
    ax = gca;
    p=imagesc(ax,domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},squeeze(plot_value(1,:,:)) );
    hold on;
    caxis(ax,[max(max(max(plot_value(:,:,:))))-dBVALUES max(max(max(plot_value(:,:,:))))] )
    if Bubble.N >5
        rectangle(ax,'Position',[domain.min_dim(domain.dimval(3)) domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(3))-domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(2))-domain.min_dim(domain.dimval(2))],'LineStyle','--','EdgeColor','white','LineWidth',2)
    elseif (Bubble.N <=5 && Bubble.N > 0)
        plot(ax,Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),'x','Color',[0.8 0.8 0.8],'MarkerSize',10,'LineWidth',3)
    end
    
    title_txt = 'Scattered Pressure Field, t_i = ';
    if (sum(sum(sum(plot_value==Linear_Field)))==numel(plot_value)) ;title_txt = 'Linear Pressure Field, t_i = ';end
    if (sum(sum(sum(plot_value==Total_Field)))==numel(plot_value)) ;title_txt = 'Total Pressure Field, t_i = ';end
    title(ax,[title_txt,num2str(1)])
    xlabel(ax,dslice.xlabel)
    ylabel(ax,dslice.ylabel)
    
    xlim(ax,[min(domain.par{domain.dimval(3)}) max(domain.par{domain.dimval(3)})])
    ylim(ax,[min(domain.par{domain.dimval(2)}) max(domain.par{domain.dimval(2)})])
    
    BubbleCluster_Colormaps(file)
    c = colorbar;
    c.Label.String = 'Pressure [dB re 1 Pa]';
    set(ax,'YDir','normal')
    if dslice.savedim(dslice.num) =='z'; set(ax,'XDir','reverse');set(ax,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
    axis image;
    set(gcf, 'Color', 'w');
    set(gca,'FontSize',30)
    
    F(1) = getframe(gcf);
    writeVideo(mov,F(1));
    
    %     hold off;
    % ======================================== Create Frames and write each frame to video =============================
    for i=2:i_end  % Change this to i_end for full time simulation movie
        set(ax.Title,'String',[title_txt,num2str(i)])
        set(p,'XData',domain.par{domain.dimval(3)},'YData',domain.par{domain.dimval(2)},'CData',squeeze(plot_value(i,:,:)))
        drawnow
        if (strcmp(file.save_movies,'yes'))
            F(i) = getframe(gcf);
            writeVideo(mov,F(i));
        end
        %         pause(0.01)
        
        
    end
    close(mov);
end