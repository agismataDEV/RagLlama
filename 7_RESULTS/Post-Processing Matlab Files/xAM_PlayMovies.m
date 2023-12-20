% CREATING MATRICES FOR NONLINEAR AND CONTRAST DATA AND PLOTTING THE PRESSURE FIELD
%
% xAM_PlayMovies(medium,domain,dslice,file,plin,pnl,Bubble)
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

function xAM_PlayMovies(medium,domain,dslice,file,plin,pnl,ScPop)

if (strcmp(file.play_movies,'yes'))
    mov = VideoWriter('GG_4Media_25kPa_1MHz_Scatter.avi','Motion JPEG AVI');
    mov.FrameRate = 20; mov.Quality = 75;
    open(mov);
    
    disp(' ')
    disp(['Play movie ...'])
    
    sizeP = size(pnl.xwave_data);
    %     Linear_Field=-500*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    %     X_Wave=-500*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    %     Right_Wave=-1*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    %     Left_Wave=-1E30*ones(sizeP(1)+sizeP(3),sizeP(2),sizeP(3));
    %%
    if (domain.dimval(1) ==1 || domain.dimval(1)==2)   % For X and Y , the domain should be moved in order to get a nice movie
        %%
        %         plin.omega       = 2*pi * 1/domain.dtpar * ((0:domain.tdimpar-1)-ceil((domain.tdimpar-1)/2))/domain.tdimpar;
        %         plin.mult = zeros(size(plin.data));  % T,Z so to create the change in z direction
        %         for i_z = 1:size(plin.data,3)
        %             plin.mult(:,:,i_z) = squeeze(repmat(exp(-1j * plin.omega .* (domain.par{3}(i_z)-domain.dxpar)*1E-3/medium.c0),1,1,size(plin.data,2)));
        %         end
        %
        %         Linear_Field = real(ifft( fft(plin.data) .* plin.mult));
        %         Linear_Field = 20*log10(abs(Linear_Field));
        for k=1:sizeP(3)
            Linear_Field(k:sizeP(1)+k-1,:,k)=20*log10(abs(pnl.xwave_data(:,:,k)));
            %             X_Wave(k:sizeP(1)+k-1,:,k)=20*log10(abs(pnl.xwave_data(:,:,k)));
            %             Right_Wave(k:sizeP(1)+k-1,:,k)=abs((pnl.rightwave_data(:,:,k)));
            %             Left_Wave(k:sizeP(1)+k-1,:,k)=abs(pnl.leftwave_data(:,:,k));
        end
        i_end = sizeP(1)+sizeP(3); % For X And Y the plot is a moving window at each Z position
        
    else   % For Z , the movie can be ploted straigh away
        for i_t=1:sizeP(1)
            Linear_Field(i_t+sizeP(3),:,:) =20*log10(abs(plin.data(i_t,:,:)));
            X_Wave(i_t+sizeP(3),:,:)=20*log10(abs(pnl.xwave_data(i_t,:,:)));
            Right_Wave(i_t+sizeP(3),:,:)=20*log10(abs(pnl.rightwave_data(i_t,:,:)));
            Left_Wave(i_t+sizeP(3),:,:) =20*log10(abs(pnl.leftwave_data(i_t,:,:)));
        end
        i_end = sizeP(1)+sizeP(3);   % For Z the plot is only for the time points
    end
    %%
    rotate = 1;
    i_end = sizeP(1) + sizeP(3);
    plot_value = squeeze(20*log10(abs(pnl.xwave_data + plin.data )));
    dBVALUES = 80;0.3*max(max(max(plot_value*1E3)));
    
    if(domain.a_t==1); plot_value = Linear_Field; end
    
    X = domain.par{domain.dimval(3)}; Y = domain.par{domain.dimval(2)};
    if (rotate);X = domain.par{domain.dimval(2)}; Y = domain.par{domain.dimval(3)} ; dimorder=[2 3 1]; A=rot90(permute(plot_value,dimorder),3);plot_value=ipermute(A,dimorder);end
    
    
    %% ========================= Initialize first frame=========================
    % This is done to speedup the process of making images , the most time is lost in making the plots
    num_frames = i_end;
    F(num_frames) = struct('cdata',[],'colormap',[]);
    
    f2 = figure('WindowState','maximized');
    ax = gca;
    p=imagesc(ax,X,Y,squeeze(plot_value(1,:,:)));
    ylabel(ax,dslice.ylabel,'Interpreter','latex')
    xlabel(ax,dslice.xlabel,'Interpreter','latex')
    if (rotate); ylabel(ax,dslice.xlabel,'Interpreter','latex'); xlabel(ax,dslice.ylabel,'Interpreter','latex'); end
    hold on;
    
%     for i=1:4
%         
%         rect_w =  (0.01*(domain.max_dim(domain.dimval(1))+domain.dypar - (domain.min_dim(domain.dimval(1))-domain.dypar)));
%         %     r = rectangle('Position',[ScPop{i}.LocRange(1,1)  ScPop{i}.LocRange(3,1) diff(ScPop{i}.LocRange(1,:)) diff(ScPop{i}.LocRange(3,:))],...
%         %     'EdgeColor',[0.64,0.08,0.18], 'FaceColor' , 'none','LineWidth',4,'LineStyle','--');%'#8BE5D3'); %% Visualization of slice
%         r = rectangle('Position',[ScPop{i}.LocRange(3,1)  ScPop{i}.LocRange(1,1) diff(ScPop{i}.LocRange(3,:))-rect_w diff(ScPop{i}.LocRange(1,:))-rect_w],...
%             'EdgeColor',[128,128,128]./255*mod(i,2)+abs(mod(i,2)-1)*[0 128 255]./255, 'FaceColor' , 'none','LineWidth',4,'LineStyle','--');%'#8BE5D3'); %% Visualization of slice
%     end
    % First MB Cloud
    Bubble=ScPop{1};
%     [xunit1,yunit1] = circle([sum(Bubble.LocRange(1,:))/2, sum(Bubble.LocRange(3,:))/2], diff(Bubble.LocRange(1,:))/2);
%     plot(xunit1,yunit1, '--','Color',"#0072BD",'LineWidth',3);
    
%     % Second MB Cloud
%     [xunit2,yunit2] = circle([sum(Bubble2.LocRange(3,:))/2, sum(Bubble2.LocRange(1,:))/2], diff(Bubble2.LocRange(3,:))/2);
%     plot(xunit2,yunit2, '--','Color',"#0072BD",'LineWidth',3);
    
    if ScPop{2}.N>0
        % LS Cloud
        Lin=ScPop{2};
        r2 = rectangle('Position',[Lin.LocRange(1,1)  Lin.LocRange(3,1) diff(Lin.LocRange(1,:)) diff(Lin.LocRange(3,:))],...
            'EdgeColor',[0.64,0.08,0.18], 'FaceColor' , 'none','LineWidth',4,'LineStyle','--');%'#8BE5D3'); %% Visualization of slice
    end
    caxis(ax,[max(max(max(plot_value)))-dBVALUES max(max(max(plot_value)))])
    title_txt = 'Scattered Pressure Field, $t_i$ = ';
    %     if (sum(sum(sum(plot_value==Linear_Field)))==numel(plot_value)) ;title_txt = 'Linear Pressure Field, t_i = ';end
    title(ax,[title_txt,num2str(0*domain.dtpar),' $\mu s$'])
    axis image;
    
    xlim(ax,[min(X) max(X)])
    ylim(ax,[min(Y) max(Y)])
    
    xAM_Colormaps(file.plot_colour)
    c = colorbar;
    c.Label.String = 'Normalized Pressure [dB]';
    c.TickLabelInterpreter='latex';
    c.Label.Interpreter='latex';
    
    %     set(ax,'YDir','normal')
    %     if dslice.savedim(dslice.num) =='z'; set(ax,'XDir','reverse');set(ax,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
    %     xlim([0 8])
    %     ylim([0 8])
    set(gcf, 'Color', 'w');
    set(gca,'FontSize',25)
    ax.TickLabelInterpreter='latex';
    %     set(gca,'LooseInset', get(gca,'TightInset'), [ 0.2 0 0 0])
    
    %     F(1) = getframe(gca);
    %     writeVideo(mov,F(1));
    %%
    % ======================================== Create Frames and write each frame to video =============================
    
    for j=2:i_end
        set(ax.Title,'String',[title_txt,num2str(round((j-1)*domain.dtpar*1E6,1)),' $\mu s$'])
        set(p,'CData',squeeze(plot_value(j,:,:)))
        drawnow
        if (strcmp(file.save_movies,'yes'))
            F(j) = getframe(gcf);
            writeVideo(mov,F(j));
        end
    end
    % j=0;
    %     for i=2:i_end  % Change this to i_end for full time simulation movie
    %         set(ax.Title,'String',[title_txt,num2str(round((i+j)*domain.dtpar*1E6)), '$\mu$s'],'Interpreter','latex')
    % %         set(p,'CData',rot90(squeeze(plot_value(i,:,:)),3))
    %         set(p,'CData',squeeze(plot_value(i,:,:)))
    % %         ylim([-1 1])
    %         drawnow
    %         if (strcmp(file.save_movies,'yes'))
    %             F(i+j) = getframe(gcf);
    %             writeVideo(mov,F(i+j));
    %         end
    %     end
    close(mov);
    
    %% =============================================== Compute angle of plane-wave =======================================
    for i=500:500 % Campure a specific time frame and plot it in a figure
        set(ax.Title,'String',[title_txt,num2str(round((i-1)*domain.dtpar*1E6,1)),' $\mu s$'])
        set(p,'CData',squeeze(plot_value(i,:,:)))
        drawnow
        drawnow
        if (strcmp(file.save_movies,'yes'))
            F(i) = getframe(gcf);
            writeVideo(mov,F(i));
        end
    end
    %%
    plot([domain.par{3}(1) domain.par{3}(end)],[0 0],'--','Color','white','LineWidth',2)
    %Capture the indexes of the max values in this plot
    [max_val,max_idx] = max(squeeze(plot_value(i,:,:))');
    
    % Plot the values in a figure
    f3 = figure('WindowState','maximized');
    plot(domain.par{3}(max_idx),domain.par{1},'LineWidth',2)
    
    xlabel(dslice.xlabel)
    ylabel(dslice.ylabel)
    axis image;
    xlim([domain.par{3}(1) domain.par{3}(end)])
    
    % Add 2 points in the same figure
    hold on;
    P_ZX = [domain.par{3}(max_idx) ; domain.par{1}]';
    P1 = P_ZX(230,:);
    P2 = P_ZX(300,:);
    plot(P1(:,1),P1(:,2),'ro','LineWidth',2','MarkerSize',5)
    plot(P2(:,1),P2(:,2),'ro','LineWidth',2','MarkerSize',5)
    
    %Find many differences
    P1_many = P_ZX(233:7:260,:);
    P2_many = P_ZX(273:7:300,:);
    Angle_many=90-atand((P2_many(:,2)-P1_many(:,2))./(P2_many(:,1)-P1_many(:,1)));
    title(['Max pressure location in X-Z plane, ' num2str(domain.angle) ' deg, @ Y = 0 mm '])
    %     title(['Average angle between the 2 points, ',num2str(mean(Angle_many)) ,' deg'])
    legend(['Angle between 2 points: ',num2str(mean(Angle_many)) ,' deg'])
    
    set(gcf, 'Color', 'w');
    set(gca,'FontSize',18)
end


% function [xunit,yunit] = circle(C,R)
% x=C(1);y=C(2);
% hold on
% th = 0:pi/50:2*pi;
% xunit = R * cos(th) + x;
% yunit = R * sin(th) + y;
% end