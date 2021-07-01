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

function [domain] = BubbleCluster_LocPlot_Paper(domain,dslice,file,Bubble)

%%
disp('Ploting the microbubble cluster...')
%=============================== 3D Domain View ==============================
f4 = figure('Position',[31,30,1288,882]);
% subplot(1,3,1)
hold on;
box on;
% if Bubble.N==1
% [X,Y,Z] = sphere;
% factor = 0.2 ; X = X*factor;Y = Y*factor;Z = Z*factor;
% s3 = surf(X+Bubble.LocGlob(:,3),Y+Bubble.LocGlob(:,1),Z+Bubble.LocGlob(:,2) , 'FaceColor','white');
% rotate(s3,[1 0 0],50)
% s3.EdgeColor ='none';
% % shading interp
% s3.FaceColor = '#404040'; s3.FaceAlpha=0.8;
% endif
[Z,X] = meshgrid(domain.par{3},domain.par{1});
s=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],dslice.pos(dslice.num)*ones(4,4) ,...
    'EdgeColor','black','FaceColor','#0080FF', 'FaceAlpha' ,0.6);  %

% text_datalabels = [ '\color{black}X \color[rgb]{0.0784 0.5725 0.7882}\bf',num2str(Bubble.LocGlob(:,1)),...
%     '\newline\color{black}Y \color[rgb]{0.0784 0.5725 0.7882}\bf',num2str(Bubble.LocGlob(:,2)),...
%     '\newline\color{black}Z \color[rgb]{0.0784 0.5725 0.7882}\bf',num2str(Bubble.LocGlob(:,3))];
% text(Bubble.LocGlob(:,3)-0.02,Bubble.LocGlob(:,1)-0.5,Bubble.LocGlob(:,2)+1.2,text_datalabels,'BackgroundColor', 'white','FontSize',20,'EdgeColor','black')
% Create textarrow
annotation('arrow',[0.12888198757764 0.428571428571429],...
    [0.780179138321995 0.886621315192744],'LineWidth',1.2)
text('Position',[-4.288198757764 0.428571428571429 4.680179138321995 ],...
    'String','Direction of Propagation','Rotation',14,'FontSize',20);

shading interp;
if (domain.dimval(1)==1)
rotate(s,[1 0 0 ], 90)
elseif (domain.dimval(1)==2)
rotate(s,[0 0 1 ], 0)
elseif (domain.dimval(1)==3)
rotate(s,[0 1 0], -90)
end
p=plot3(Bubble.LocGlob(:,3),Bubble.LocGlob(:,1),Bubble.LocGlob(:,2),'o','MarkerSize',7/(log10(Bubble.N)+1),'MarkerFaceColor','black');

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

s=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],dslice.pos(dslice.num)*ones(4,4) ,...
    'EdgeColor','black','FaceColor','#0080FF', 'FaceAlpha' ,0.6);  %

if (strcmp(file.plot_colour,'gray')); s.FaceColor =[ 192 192 192]/255;s.EdgeColor='black'; p.MarkerEdgeColor= '#404040' ;end
% 
% s1=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],Bubble.LocGlob(2)*ones(4,4) ,...
%     'EdgeColor','#404040','FaceColor','#36EDFA', 'FaceAlpha' ,0.3);
xlim([min(domain.par{3}) max(domain.par{3})])
ylim([min(domain.par{1}) max(domain.par{1})])
zlim([min(domain.par{2}) max(domain.par{2})])

x = xlabel('Z [mm]');
set(x, 'Position',get(x, 'Position').*[1,1,1],'Rotation',10)
y = ylabel('X [mm]');
set(y, 'Position',get(y, 'Position').*[1,1,1],'Rotation',-21)
zlabel('Y [mm]')
title(['Computational Domain with ',num2str(Bubble.N),' Scatterer(s) /ml'])
grid on
grid minor
hold off;
set(gca,'FontSize',25)
set(gcf,'Color','white')
%% ======================= 3D Microbubble Cluster Focused position =============
% subplot(1,2,1)
% f4 = figure('WindowState','maximized');

f5=figure;
set(f5,'Position',2*get(f5,'Position'));
hold on;
box on;
[X,Y,Z] = meshgrid(domain.par{domain.dimval(2)},domain.par{domain.dimval(1)},domain.par{domain.dimval(3)});
% s=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],dslice.pos(dslice.num)*ones(4,4),...
%     'EdgeColor','black', 'FaceAlpha' ,0.6);
zlim_foc = (round(10*([min(Bubble.LocGlob(:,domain.dimval(3))) max(Bubble.LocGlob(:,domain.dimval(3)))] + [-domain.dzpar domain.dzpar]))/10).*[1 1];
xlim_foc = round(10*([min(Bubble.LocGlob(:,domain.dimval(2))) max(Bubble.LocGlob(:,domain.dimval(2)))] + [-domain.dxpar domain.dxpar]))/10;
ylim_foc = (round(10*([min(Bubble.LocGlob(:,domain.dimval(1))) max(Bubble.LocGlob(:,domain.dimval(1)))] + [-domain.dypar domain.dypar]))/10).*[ 1 1];

p = plot3(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(2)),Bubble.LocGlob(:,domain.dimval(1)),'o','MarkerSize',8/(log10(Bubble.N)+1));
s = plot3([zlim_foc(1) zlim_foc(1) zlim_foc(end) zlim_foc(end) zlim_foc(1)],xlim_foc(1)*ones(5,1), [ylim_foc(1) ylim_foc(end) ylim_foc(end) ylim_foc(1) ylim_foc(1)], 'Color','black','LineWidth',3);
if (strcmp(file.plot_colour,'gray')); p.MarkerEdgeColor= 'black';p.MarkerFaceColor='black';end%'#404040' ;end
view(-37.5,30)

zlim_foc3D = (round(10*([min(Bubble.LocGlob(:,domain.dimval(3))) max(Bubble.LocGlob(:,domain.dimval(3)))] + [-domain.dzpar domain.dzpar]))/10);
xlim_foc3D = round(10*([min(Bubble.LocGlob(:,domain.dimval(2))) max(Bubble.LocGlob(:,domain.dimval(2)))] + [-domain.dxpar domain.dxpar]))/10;
ylim_foc3D = (round(10*([min(Bubble.LocGlob(:,domain.dimval(1))) max(Bubble.LocGlob(:,domain.dimval(1)))] + [-domain.dypar domain.dypar]))/10);
xlim(zlim_foc3D)
ylim(xlim_foc3D)
zlim(ylim_foc3D)

x = xlabel(dslice.xlabel);
set(x, 'Position',get(x, 'Position').*[1,1,1.1],'Rotation',10)
y = ylabel(dslice.ylabel);
set(y, 'Position',get(y, 'Position').*[1,1,1.1],'Rotation',-21)
zlabel(dslice.zlabel)
title('3D View of Cloud')
grid on
hold off;

set(gca,'FontSize',30)
set(gcf, 'Color', 'w');
%====================== X-Z Side View of Focused Region ==========================
% subplot(1,2,2)
% f5 = figure('WindowState','maximized');
f6=figure;
set(f6,'Position',2*get(f6,'Position'));
hold on;
p=plot(Bubble.LocGlob(:,domain.dimval(3)),Bubble.LocGlob(:,domain.dimval(1)),'o','MarkerSize',8/(log10(Bubble.N)+1));

% Create dimensions that encloses the microbubble cluster but includes the
% slice of the domain stored in the output file
domain.min_dim = min([Bubble.LocGlob ; domain.par{1}(end) domain.par{2}(end) domain.par{3}(end)]);
domain.max_dim = max([Bubble.LocGlob ; domain.par{1}(1) domain.par{2}(1) domain.par{3}(1)]);
domain.min_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)) = min(domain.min_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)), dslice.pos(dslice.num));
domain.max_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)) = max(domain.max_dim(['x' 'y' 'z']==dslice.savedim(dslice.num)), dslice.pos(dslice.num));

rect_w =  0.01*(domain.max_dim(domain.dimval(1))+domain.dypar - (domain.min_dim(domain.dimval(1))-domain.dypar));
% r = rectangle('Position',[floor(10*(domain.min_dim(domain.dimval(3))-domain.dzpar))/10 dslice.pos(dslice.num)-rect_w  ceil(10*(domain.max_dim(domain.dimval(3))+domain.dzpar-(domain.min_dim(domain.dimval(3))-domain.dzpar)))/10 rect_w],...
%     'EdgeColor','k', 'FaceColor' , '#8BE5D3'); %% Visualization of slice
if (strcmp(file.plot_colour,'gray')); p.MarkerEdgeColor = 'black';p.MarkerFaceColor='black';end%'#404040' ;end%; r.FaceColor = [ 192 192 192]/255;end

% xlim_out = round(10*([ domain.min_dim(domain.dimval(3)) domain.max_dim(domain.dimval(3))] + [-domain.dzpar domain.dzpar]))/10; 
% ylim_out = round(10*([ domain.min_dim(domain.dimval(1)) domain.max_dim(domain.dimval(1))] + [-domain.dypar domain.dypar]))/10;
xlim_out = zlim_foc;
ylim_out = ylim_foc;
% xlim(sum(xlim_out)/2 + pdist(xlim_out')/1*[-1 1])
% ylim(sum(ylim_out)/2 + pdist(ylim_out')/1*[-1 1])
xlim(xlim_out)
ylim(ylim_out)

xlabel(dslice.xlabel)
ylabel(dslice.zlabel)

title([dslice.zlabel(1),' - ',dslice.xlabel(1),' Plane View of Cloud'])
grid on
hold off;
set(gca,'FontSize',30)
set(gcf, 'Color', 'w');
%============================= Create subplot title ======================
if (strcmp(file.saveplot,'yes'))
        figure(f5)
%     saveas(f4,[file.savedir,'/',file.dirname,'_3D_cluster_slice',int2str(dslice.num),'.jpg'])
        folder = ['../',file.savedir,'/',file.dirname,'_3D_cluster_slice',int2str(dslice.num),'.png'];
        cd export_fig/
        export_fig  folder -painters -q110 
        movefile('folder.png',folder)
        cd ../
%     saveas(f5,[file.savedir,'/',file.dirname,'_YZ_cluster_slice',int2str(dslice.num),'.jpg'])
        figure(f6)
    folder1 = ['../',file.savedir,'/',file.dirname,'_YZ_cluster_slice',int2str(dslice.num),'.png'];
        cd export_fig/
        export_fig  folder1 -painters -q110
        movefile('folder1.png',folder1)
        cd ../
end
%%

%=============================== 3D Domain View ==============================
f4 = figure('Position',[31,30,1288,882]);
LocPlot_FontSize = 20;
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
p=plot3(Bubble.LocGlob(:,3),Bubble.LocGlob(:,1),Bubble.LocGlob(:,2),'o','MarkerSize',7/(log10(Bubble.N)+1),'MarkerFaceColor','black');
if (strcmp(file.plot_colour,'gray')); s.FaceColor =[ 192 192 192]/255;s.EdgeColor='black'; p.MarkerEdgeColor= '#404040' ;end

annotation('arrow',[0.12888198757764 0.428571428571429],...
    [0.780179138321995 0.886621315192744],'LineWidth',1.2)
text('Position',[-8.344848181658563,7.004384680729004,7.244327405467004],...
    'String','Direction of Propagation','Rotation',14,'FontSize',20);

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
% title(['Microbubble Cluster of ',integerPart,' Bubble(s) per ',num2str(round(Bubble.Volume*1e-3,2)),' [mL] at ',num2str(dslice.savedim(dslice.num)), ' = ', num2str(dslice.pos(dslice.num)) , ' [mm]'])
grid on
hold off;

set(gca,'FontSize',LocPlot_FontSize)
set(gcf, 'Color', 'white');
%======================= 3D Microbubble Cluster Focused position =============
f5= figure('Position',[31,30,1288,882]);
hold on;
box on;
[Z,X] = meshgrid(domain.par{domain.dimval(3)},domain.par{domain.dimval(2)});
s=surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],dslice.pos(dslice.num)*ones(4,4),...
    'EdgeColor','black', 'FaceAlpha' ,0.6);
p=plot3(Bubble.LocGlob(:,3),Bubble.LocGlob(:,1),Bubble.LocGlob(:,2),'o','MarkerSize',7/(log10(Bubble.N)+1),'MarkerFaceColor','black');
if (strcmp(file.plot_colour,'gray')); s.FaceColor =[ 192 192 192]/255;s.EdgeColor='black'; p.MarkerEdgeColor= '#404040' ;end
view(-37.5,30)

%
xlim(round(10*([min(Bubble.LocGlob(:,domain.dimval(3))) max(Bubble.LocGlob(:,domain.dimval(3)))] + [-domain.dzpar domain.dzpar]))/10)
ylim(round(10*([min(Bubble.LocGlob(:,domain.dimval(2))) max(Bubble.LocGlob(:,domain.dimval(2)))] + [-domain.dxpar domain.dxpar]))/10)
zlim(round(10*([min(Bubble.LocGlob(:,domain.dimval(1))) max(Bubble.LocGlob(:,domain.dimval(1)))] + [-domain.dypar domain.dypar]))/10)

x = xlabel(dslice.xlabel);
set(x, 'Position',get(x, 'Position').*[1,1,1],'Rotation',10)
y = ylabel(dslice.ylabel);
set(y, 'Position',get(y, 'Position').*[1,1,1],'Rotation',-21)
zlabel(dslice.zlabel)
title('Microbubble Cluster Focused Region')
grid on

annotation('arrow',[0.12888198757764 0.428571428571429], [0.780179138321995 0.886621315192744],'LineWidth',1.2)
text('Position',[6.918560488193179,-2.760877433698681,4.718873489954944 ],'String','Direction of Propagation','Rotation',14,'FontSize',20);
hold off;

set(gca,'FontSize',LocPlot_FontSize)
set(gcf, 'Color', 'white');
%========================== Save plot ===================================
if (strcmp(file.saveplot,'yes'))
        figure(f4)
%     saveas(f4,[file.savedir,'/',file.dirname,'_3D_cluster_slice',int2str(dslice.num),'.jpg'])
        folder = ['../',file.savedir,'/',file.dirname,'_3D_cluster_slice',int2str(dslice.num),'.png'];
        cd export_fig/
        export_fig  folder -painters -q110 
        movefile('folder.png',folder)
        cd ../
%     saveas(f5,[file.savedir,'/',file.dirname,'_YZ_cluster_slice',int2str(dslice.num),'.jpg'])
        figure(f5)
    folder1 = ['../',file.savedir,'/',file.dirname,'_YZ_cluster_slice',int2str(dslice.num),'.png'];
        cd export_fig/
        export_fig  folder1 -painters -q110
        movefile('folder1.png',folder1)
        cd ../
end
end