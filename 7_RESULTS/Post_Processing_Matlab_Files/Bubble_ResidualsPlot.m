% ***********************************************************************************************************
% PLOT THE CONVERGENCE ERROR OF THE CODE BASED ON NUMBER OF SCATTERERS
%
% *************************************************************************
% Bubble_ResidualsPlot(medium,domain,slice,file,plin,pnl,Bubble,BubbleList)
% *************************************************************************
% Plot the convergence error for different modes of cluster by changing the BubbleList
% input which describes the number of the scatterers
%
% Inputs:
% medium        struct      Contains all the medium properties, e.g Density, Speed of Sound etc
% domain        struct      Contains all the computational domain properties, e.g space and time step
% slice         struct      Contains all the information about the slice ,e.g Dimension, Position, Index etc
% file          struct      Contains all the loaded file properties ,e.g directory and root name etc
% plin          struct      Contains all the pressure values of the plin slice
% pnl           struct      Contains all the pressure values of the pnl slice
% Bubbles       struct      Contains all properties of the Scatterers' Cluster
% BubbleList    integer     Array that contains different number for Bubble Cluster
%                           In this way, the code can compare convergence error for different modes
%
% ******************
% AM, version 200723
% ******************
%
% ***********************************************************************************************************


function Bubble_ResidualsPlot(medium,domain,dslice,file,plin,pnl,Bubble,BubbleList,filename)

disp(['Ploting the Residuals ... '])


% % plots
% figure(1)
% Bubble.N = [1,10,100,1000, 10000,100000];
% % hold on
% for Bubble.N = 1:length(Bubble.N)
% name=[int2str(Bubble.N(Bubble.N)),'_qsource',concentrated,'\Neumann'] ;                  % file name
% % name='BiCGSTAB'                  % file name
% % name='SteepestDescent'           % file name
%
% ProcessorsNumber=8;              % number of processors used for the simulation
%
% for i=0:ProcessorsNumber-1
%     if i<10
%     err = textread(['.\',name,'errorN000000000',num2str(i),'.dat']);
%     else
%     err = textread(['.\',name,'errorN00000000',num2str(i),'.dat']);
%     end
%     Dim=size(err);
%     TotalError(:,i+1)=err(:,1);
% end
%
% error=sum(TotalError');
%
% for i=1:Dim(1)-1
%     Residual(i+1)=sqrt(error(i+1)./error(1));
% end
% Residual(1)=1;
%
% i_end = 's';
% if Bubble.N==1 ; i_end='';end
% legend_data1{Bubble.N} = [int2str(Bubble.N(Bubble.N)), ' Bubble',i_end ' /1st Iter'];
% Error1(Bubble.N,:) = Residual;
% end
%
%
% % Create multiple lines using matrix input to plot
% plot1 = semilogy(0:Dim(1)-1,Error1);
% % axes1 = axes;
% grid on;
% set(plot1(1),'DisplayName','1 Bubble','Marker','o','LineWidth',1,...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
% set(plot1(2),'DisplayName','10 Bubbles','Marker','x',...
%      'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
% set(plot1(3),'DisplayName','100 Bubbles','Marker','square',...
%      'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
% set(plot1(4),'DisplayName','1000 Bubbles','Marker','diamond',...
%      'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
%
% % Create ylabel
% ylabel('residual log[Err(n)]');
%
% % Create xlabel
% xlabel('Number of Iterations [n]');
% % Set the remaining axes properties
% % legend(legend_data1)

%% ERR1(n) analisys
% Create axes

if (strcmp(file.plot_converr,'yes'))
    for Bubblei = 1:1
        LS = split(cellstr(ls([file.dirname,'/',file.rootname,'*_',dslice.savedim(dslice.num),int2string_ICS(dslice.num),int2string_ICS(0),int2string_ICS(0),'*.h5'])),[file.rootname,'_']);
        
        LS = split(LS(:,2) ,'_y_');
        IterationNumber = max(str2double(LS(:,1))); IterationNumber=25;
        
        ii = 0;
        txt_err = {'total'; 'contrast'};
%         file_dirname = [filename, erase(num2str(BubbleList(Bubblei),'%.0E'),'+0'),'MBs'];
        file_dirname = file.dirname;
        
        for i =1:1
            if (i==1)
                current   = [file_dirname '/' file.rootname int2string_ICS(0) '_' file.focalplanename int2string_ICS(ii)];
            else
                current   = [file_dirname '/' file.contrast_name int2string_ICS(1) '_' file.focalplanename int2string_ICS(ii)];
            end
            output = load_ICS_slice_par(current);
            p_curr  = squeeze(output.data); p_first = p_curr;
                
            for jj=18:IterationNumber+mod(i,2)
                    p_prev = p_curr; 
                    prev = current;
                    disp(['Loading ... ',current])
                    
                    filename2 = file.rootname; if i==2;filename2 = file.contrast_name;end
                    current       = [file_dirname '/' filename2 int2string_ICS(jj-mod(i,2)) '_' file.focalplanename int2string_ICS(ii)];
                    output = load_ICS_slice_par(current);
                    p_curr  = squeeze(output.data);

                    if (jj==2) ;eval(['dp = p_curr-p_prev;']);end
                    Error_prev{i}(Bubblei,jj-1) = sqrt(sum(sum(sum((p_curr-p_prev).^2)))./sum(sum(sum(p_prev.^2))));
                    Error_dp{i}(Bubblei,jj-1) =  sqrt(sum(sum(sum((p_curr-p_prev).^2)))./sum(sum(sum(dp.^2))));
                    Error_init{i}(Bubblei,jj-1) = sqrt(sum(sum(sum((p_curr-p_prev).^2)))./sum(sum(sum(p_first.^2))));
                    dp_val{i}(Bubblei,jj-1) = squeeze(max(abs(p_curr(:,round(size(p_curr,2)/2)+1,470)-p_prev(:,round(size(p_curr,2)/2)+1,470))));
%                     eval(['Error_prev_',txt_err{i},'(Bubblei,jj-1) = sqrt(sum(sum(sum((p_curr-p_prev).^2)))./sum(sum(sum(p_prev.^2))));']);
%                     eval(['Error_dp_',txt_err{i},'(Bubblei,jj-1) = sqrt(sum(sum(sum((p_curr-p_prev).^2)))./sum(sum(sum(dp.^2))));']);
%                     eval(['Error_init_',txt_err{i},'(Bubblei,jj-1) = sqrt(sum(sum(sum((p_curr-p_prev).^2)))./sum(sum(sum(p_first.^2))));']);
%                     eval(['dp_',txt_err{i},'(Bubblei,jj-1) = squeeze(max(abs(p_curr(:,round(size(p_curr,2)/2)+1,470)-p_prev(:,round(size(p_curr,2)/2)+1,470))));']);
            end
            Error_prev{i}(Bubblei,1:jj) = ([1,Error_prev(Bubblei,2:jj)]);
            Error_init{i}(Bubblei,1:jj) = ([1,Error_init(Bubblei,2:jj)]);
            Error_dp{i}(Bubblei,1:jj) = ([1,Error_dp(Bubblei,2:jj)]);
        end
        
        % PLOTS
        i_end = 's';
        legend_data2{Bubblei,1} = [int2str(BubbleList(Bubblei)), ' MB',i_end ];
        legend_data2{Bubblei,2} = [int2str(BubbleList(Bubblei)), ' MB',i_end ];
        %  plot(0:IterationNumber,log([1,Error]))
    end
    %% Plot
    % hold on;
    % Create multiple lines using matrix input to plot
%     hold on
        f8=figure('WindowState','maximized');
    for i = 1:1
        hold on;
        for Bubblei =1:4
               plot_h(Bubblei) =  semilogy(0:size(Error_prev{i}(Bubblei,:),2)-1,flipud(dp_val{i}(Bubblei,:)));
    %         eval(['plot',num2str(i),' = semilogy(0:size(Error_prev_',txt_err{i},',2)-1,[flipud(dp_',txt_err{i},')])'])
    %            
    %         set(eval(['plot',num2str(i),'(1)']),'DisplayName',[num2str(BubbleList(1)),' Bubble(s) Pres'],'Marker','o','MarkerSize',12,'LineWidth',2,'LineStyle','--',...
    %             'Color','#994C00');
    %         set(eval(['plot',num2str(i),'(2)']),'DisplayName',[num2str(BubbleList(2)),' Bubble(s) Pres'],'Marker','x','MarkerSize',12,'LineWidth',2,'LineStyle','--',...
    %             'Color','#606060');
    %         set(eval(['plot',num2str(i),'(3)']),'DisplayName',[num2str(BubbleList(3)),' Bubble(s) Pres'],'Marker','*','MarkerSize',12,'LineWidth',2,'LineStyle','--',...
    %             'Color','#009999');
    %         set(eval(['plot',num2str(i),'(4)']),'DisplayName',[num2str(BubbleList(4)),' Bubble(s) Pres'],'Marker','square','MarkerSize',12,'LineWidth',2,'LineStyle','--',...
    %             'Color','#990000');
    %         set(plot2(4),'DisplayName','1000 Bubbles','Marker','diamond',...
    %             'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
    %         set(plot2(5),'DisplayName','10000 Bubbles','Marker','o',...
    %             'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
    %         set(plot2(6),'DisplayName','100000 Bubbles','Marker','*',...
    %             'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
        end
        grid on;
        % Create ylabel
        ylabel('RRMSE(n)');
        % Create xlabel
        xlabel('Number of Iterations [n]');
        title(['Convergence error of ',lower(txt_err{i}),' pressure'])
        % grid(axes1,'on');
        % Set the remaining axes properties
%         legend(legend_data2{:,1})
        % legend([legend_data1 legend_data2])
        set(gca,'FontSize',30)
        set(gcf,'Color','w')    
        
        if (strcmp(file.saveplot,'yes'))
            %             saveas(f8,[file.savedir,'/',file.dirname,'_Convergence_Error_',txt_err{i},'.jpg'])
            folder1 = ['../',file.savedir,'/',file.dirname,'_Convergence_Error_',txt_err{i},'.jpg'];
            cd export_fig/
            export_fig  img1 -painters -q110
            movefile('img1.png', folder1)
            cd ../
        end
    end
%     hold off

end

if strcmp(file.plot_attenslices, 'yes')
    f9 = figure('WindowState','maximized');
    grid on;
    hold on;
    [Y,X] = meshgrid(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)});
    Z = zeros(size(X)) + dslice.pos(dslice.num);
    surf_atten = 20*log10(pnl.atten_fund./plin.atten_fund);
    H = surf(Z,X,Y,flipud(surf_atten));
    set(H,'edgecolor','none')
    colormap(flipud(parula))
    
    plot3( dslice.pos(dslice.num)*ones(1,5), ...
        [domain.min_dim(domain.dimval(2)) domain.min_dim(domain.dimval(2)) domain.max_dim(domain.dimval(2)) domain.max_dim(domain.dimval(2)) domain.min_dim(domain.dimval(2))], ...
        [domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(3)) domain.max_dim(domain.dimval(3)) domain.min_dim(domain.dimval(3)) domain.min_dim(domain.dimval(3))], ...
        '--','LineWidth',2,'Color', 'white')
    axis image
    xlim([Bubble.LocRange(3,1)-5 Bubble.LocRange(3,2)])
    zlim([domain.par{domain.dimval(3)}(1) domain.par{domain.dimval(3)}(end)])
    ylim([domain.par{domain.dimval(2)}(1) domain.par{domain.dimval(2)}(end)])
    title(['Attenuation [dB] at ',num2str(dslice.savedim(dslice.num)), ' = ', num2str(dslice.pos(dslice.num)) , ' [mm]'])
    xlabel(dslice.zlabel)
    zlabel(dslice.xlabel)
    ylabel(dslice.ylabel)
    caxis([-10 1])
    colorbar;
    view(3)
    rotate3d
    if (strcmp(file.saveplot,'yes'))
        saveas(f9,[file.savedir,'/',file.dirname,'_attenuation_surf_slice',int2str(dslice.num),'.jpg'])
    end
end
end