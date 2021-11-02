%% A.M. 23-02-2020 - This script reads and displays the results contained in the HDF5 data
clear all
% close all
clc

%% INITIALIZATION
freq0               = 1E6;                    %fundamental frequency
c0                  = 1482;                     % Speed of Sound [m/sec]
rootname            = 'TESTNeumann';
BubbleN             = 10000;
Fnyq                = 4.5;
dirname             = [int2str(BubbleN),'_',int2str(Fnyq*10),'FNYQ_NFBU_FA'];
dirname             =  [int2str(BubbleN),'_test']
dirname0            = split(dirname,'_');
beamiterations      = [1];
slicesavedim        = ['y'];
iDimXYZ             = 75;
slicesavepos        = [3.128];                  %Unormalized with the starting position
strides             = [1 1 1 1];
offsets             = [0 0 0 0];
ii                  = 0;
slicenum            = 2;
iProcN              = 48;
BubbleN             = 1000;

focalplanename = [ slicesavedim(1) int2string_ICS(size(slicesavepos,1)) ];
% transversalplanename = [ slicesavedim(2) int2string_ICS(slicesavepos(2)) ];

%% LINEAR DATA WITHOUT ANY BUBBLES
filename0 = [dirname '/' rootname ...
    int2string_ICS(0)...
    '_' focalplanename ...
    int2string_ICS(0) ];
i_start=1;
if isempty(dirname) i_start = 2; end


%% LINEAR DATA
output      = load_ICS_slice(filename0(i_start:end),strides,offsets);
plin        = squeeze(output.data);
tdimpar     = size(plin,1);
xdimpar     = size(plin,2);
zdimpar     = size(plin,3);
ydimpar     = size(plin,3); if slicesavedim(1)  =='y' ;ydimpar = iDimXYZ;end
TXYoff      = output.TXYoff;
TXYstart    = output.start(1:3);
dtpar       = output.stepsize(1);
dxpar       = output.stepsize(2)*1e+3; % step size dx in mm
dypar       = output.stepsize(3)*1e+3; % step size dy in mm
dzpar       = output.stepsize(4)*1e+3; % step size dz in mm
zstart      = output.start(4)*dzpar;
xsource     = (output.start(2)+(0:size(plin,2)-1))*dxpar;
ysource     = (output.start(3)+(0:size(plin,3)-1))*dypar;
if slicesavedim(1) =='y' iStartY = -ceil(ydimpar/2); end

%% generate axes for focal plane
TXYstarts=TXYstart+min(TXYoff);
tdimpar_glob=tdimpar+max(TXYoff(:,1))-min(TXYoff(:,1));
xdimpar_glob=xdimpar+max(TXYoff(:,2))-min(TXYoff(:,2));
tpar=(0:tdimpar-1)*dtpar;
tpar_glob=(TXYstart(1)+min(TXYoff(:,1))+(0:tdimpar_glob-1))*dtpar;
xpar=(TXYstart(2)+TXYoff(1,2)+(0:xdimpar-1))*dxpar;
xpar_glob=(TXYstart(2)+min(TXYoff(:,2))+(0:xdimpar_glob-1))*dxpar;
zpar = (zstart + (0:zdimpar-1)*dzpar);


%% NONLINEAR DATA
third       = [dirname '/' rootname '_00' num2str(beamiterations) '_y_00',num2str(slicenum),'_00' num2str(ii)];

if (beamiterations>9)&&(beamiterations<100)
    third       = [dirname '/' rootname '_0' num2str(beamiterations) '_y_00',num2str(slicenum),'_00' num2str(ii)];
else if beamiterations>99
        third       = [dirname '/' rootname '_' num2str(beamiterations) '_y_00',num2str(slicenum),'_00' num2str(ii)];
    end
end

output      = load_ICS_slice(third);
pnonl       = squeeze(output.data);
xdimpars    = xdimpar+max(TXYoff(:,2))-min(TXYoff(:,2));
ydimpars    = ydimpar+max(TXYoff(:,3))-min(TXYoff(:,3));
ypar        = (TXYstart(3)+TXYoff(1,3)+(0:ydimpar-1))*dypar;
ypars       = (TXYstart(3)+min(TXYoff(:,3))+(0:ydimpars-1))*dypar;
tpar        = (0:tdimpar-1)*dtpar;
xpar        = (TXYstart(2)+TXYoff(1,2)+(0:xdimpar-1))*dxpar;
xpars       = (TXYstart(2)+min(TXYoff(:,2))+(0:xdimpars-1))*dxpar;
zpar        = (zstart + (0:zdimpar-1)*dxpar);
ypar        = (iStartY*dypar:dypar:(ydimpar-1+iStartY)*dypar);

Dimensions=size(pnonl);
sizeP=size(pnonl);
%
%% CREATING MATRICES
plin_fund = nan(xdimpar_glob,zdimpar);
pnl_2nd  = nan(xdimpar_glob,zdimpar);
pnl_3rd  = nan(xdimpar_glob,zdimpar);
pnl_4th  = nan(xdimpar_glob,zdimpar);
pnl_5th  = nan(xdimpar_glob,zdimpar);
pnl_SH   = nan(xdimpar_glob,zdimpar);

% generate fundamental, higher harmonic and superharmonic profiles
for i=1:zdimpar
    plin_fund((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(plin(:,:,i)),dtpar,.6*freq0,1.4*freq0));
    pnl_2nd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,1.6*freq0,2.4*freq0));
    pnl_3rd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,2.6*freq0,3.4*freq0));
   pnl_4th((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,3.6*freq0,4.4*freq0));
%         pnl_5th((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(plin(:,:,i)),dtpar,4.7*freq0,5.3*freq0));
%         pnl_SH((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(plin(:,:,i)),dtpar,3.5*freq0,4.4*freq0));
end

maxfund=max(max(plin_fund));


%% PLOTS
dimensions = size(plin_fund);
figure(1)
hold on
plot(zpar-zpar(1),plin_fund( round(dimensions(1)/2),:),'k')
plot(zpar-zpar(1),pnl_2nd( round(dimensions(1)/2),:),'b')
plot(zpar-zpar(1),pnl_3rd(  round(dimensions(1)/2),:),'r')
plot(zpar-zpar(1),pnl_4th(  round(dimensions(1)/2),:),'y')
% plot(zpar-zpar(1),pnl_5th(  round(dimensions(1)/2),:),'g')
% plot(zpar-zpar(1),pnl_SH(  round(dimensions(1)/2),:),'g')
grid on
xlabel('z [m]')
ylabel('Pressure [dB re 1 Pa] ')
% legend('Linear Field', '2H', '3H', '4H', '5H')
legend('Linear Field', '2H', '3H', '4H')
% legend('Linear Field', '2H', 'SH')
% axis([0 20 50 140 ])


Nsubplot = 4;
splot_div =2;
dBVALUES=25;
figure;
subplot(Nsubplot/splot_div,Nsubplot/(Nsubplot/splot_div),1)
imagesc(zpar,xpar,plin_fund)
caxis([max(max(plin_fund))-dBVALUES max(max(plin_fund))])
axis image
title('F0')
xlabel('z (mm)')
ylabel('x (mm)')
c = colorbar;
c.Label.String = 'Pressure [dB]';

subplot(Nsubplot/splot_div,Nsubplot/(Nsubplot/splot_div),2)
imagesc(zpar,xpar,pnl_2nd)
caxis([max(max(pnl_2nd))-dBVALUES max(max(pnl_2nd))])
title('2H')
xlabel('z (mm)')
ylabel('x (mm)')
axis image
c = colorbar;
c.Label.String = 'Pressure [dB]';

subplot(Nsubplot/splot_div,Nsubplot/(Nsubplot/splot_div),3)
imagesc(zpar,xpar,pnl_3rd)
caxis([max(max(pnl_3rd))-dBVALUES max(max(pnl_3rd))])
title('3H')
xlabel('z (mm)')
ylabel('x (mm)')
axis image
c = colorbar;
c.Label.String = 'Pressure [dB]';

subplot(Nsubplot/splot_div,Nsubplot/(Nsubplot/splot_div),4)
imagesc(zpar,xpar,pnl_4th)
caxis([max(max(pnl_4th))-dBVALUES max(max(pnl_4th))])
title('4H')
xlabel('z (mm)')
ylabel('x (mm)')
axis image
c = colorbar;
c.Label.String = 'Pressure [dB]';
% figure(2)
% subplot(harm,1,5)
% imagesc(zpar,xpar,pnl_5th)
% caxis([max(max(pnl_5th))-dBVALUES max(max(pnl_5th))])
% title('5H')
% xlabel('z (mm)')
% ylabel('x (mm)')
% axis image
% colorbar
% 
% figure(2)
% subplot(harm,1,3)
% imagesc(zpar,xpar,pnl_SH)
% caxis([max(max(pnl_SH))-dBVALUES max(max(pnl_SH))])
% title('SH')
% xlabel('z (mm)')
% ylabel('x (mm)')
% axis image
% c = colorbar;
% c.Label.String = 'Pressure [dB]';

%% MICROBUBBLES
if BubbleN ==0 plin_3rd = pnl_3rd; plin_4th = pnl_4th; return;end
dLambdaNN = c0*1e+3/freq0;
for iProcID =0:iProcN-1
    
    if iProcID<10
        filename_loc = [dirname '/Bubbles/output_file_00',int2str(iProcID),'.bin'];
        filename_con = [dirname '/Bubbles/output_file_contrast_00',int2str(iProcID),'.bin'];
    else
        filename_loc = [dirname '/Bubbles/output_file_0',int2str(iProcID),'.bin'];
        filename_con = [dirname '/Bubbles/output_file_contrast_0',int2str(iProcID),'.bin'];
    end
    fileID_loc = fopen(filename_loc,'r');
    
    formatSpec = '%f %f %f %f %f %f';
    size_loc = [3 Inf];
    BubbleLoc(1+BubbleN*iProcID:(iProcID+1)*BubbleN,1:3) = fscanf(fileID_loc,formatSpec,size_loc)';
    
    fileID_con = fopen(filename_con,'r');
    size_con = [1 Inf];
    BubbleCon(iProcID*tdimpar*BubbleN+1:(iProcID+1)*tdimpar*BubbleN) = fscanf(fileID_con,formatSpec,size_con);
    
    if iProcID >0
        if sum(sum(BubbleCon(iProcID*tdimpar*BubbleN+1:(iProcID+1)*tdimpar*BubbleN)==BubbleCon((iProcID-1)*tdimpar*BubbleN+1:((iProcID-1)+1)*tdimpar*BubbleN)))~=BubbleN*tdimpar disp('ERROR in Con'); return; end
        if sum(sum(BubbleLoc(1+BubbleN*iProcID:(iProcID+1)*BubbleN,1:3)==BubbleLoc(1+BubbleN*(iProcID-1):((iProcID-1)+1)*BubbleN,1:3)))~=BubbleN*3 disp('ERROR in Loc'); return; end
    end
end

% plin_bubble = nan(tdimpar,xdimpar,zdimpar);
BubbleCon_final = BubbleCon(1:tdimpar*BubbleN);
BubbleLoc_finalN = BubbleLoc(1:BubbleN,1:3);
BubbleLoc_final = BubbleLoc_finalN*dLambdaNN;
% plin_bubble = plin(:,find(BubbleLoc_final(1)+TXYstart(2)*dxpar<xpar,1),find(BubbleLoc_final(3)+1*dzpar<zpar,1));


%%
figure;
hold on;
[Z,X] = meshgrid(zpar,xpar);
%S = surf(Z,X,slicesavepos*ones(size(X,1),size(X,2)))
S = surf([Z(1,1) Z(1,end) Z(1,1) Z(end,1)] , [X(1,1) X(1,end) X(1,1) X(end,1)],[slicesavepos(1)*ones(4,4)]);
S.AlphaData = 0;
S.FaceAlpha = 0.3;
BubbleLoc_finalG(:,1) = BubbleLoc_final(:,1)+TXYstart(2)*dxpar;
BubbleLoc_finalG(:,2) = BubbleLoc_final(:,2)+min(ypar);
BubbleLoc_finalG(:,3) = BubbleLoc_final(:,3);
plot3(BubbleLoc_finalG(:,3),BubbleLoc_finalG(:,1),BubbleLoc_finalG(:,2),'o','MarkerSize',3)
%patch([max(zpar) min(zpar) min(zpar) max(zpar)],[max(zpar) min(zpar) min(zpar) max(zpar)],[max(ypar) min(ypar)  min(ypar) max(ypar) ],...
%     'FaceColor',[0.941176470588235 0.941176470588235 0.941176470588235]);
view(-37.5,30)
% ylim([min(BubbleLoc_finalG(:,1)) max(BubbleLoc_finalG(:,1))])
% zlim([min(BubbleLoc_finalG(:,2)) max(BubbleLoc_finalG(:,2))])
% xlim([min(BubbleLoc_finalG(:,3)) max(BubbleLoc_finalG(:,3))])
ylim([min(xpar) max(xpar)])
zlim([min(ypar) max(ypar)])
xlim([min(zpar) max(zpar)])

xlabel('z [mm]')
ylabel('x [mm]')
zlabel('y [mm]')
title('Microbubble Position')
grid on
hold off;
fclose('all')


%%


xnorm = (xpar-xpar(1))/dLambdaNN;
ynorm = (0:dypar:31*dypar)/dLambdaNN;
znorm = (zpar-zpar(1))/dLambdaNN;

[X,Y,Z] = meshgrid(xnorm,ynorm, znorm);
dRadius = sqrt((BubbleLoc_finalN(1) - X).^2+(BubbleLoc_finalN(3) - Z).^2+(BubbleLoc_finalN(2)-Y).^2 );
dRadius = (reshape(dRadius,[1 size(dRadius,1) size(dRadius,2) size(dRadius,3)]));
i = floor(size(BubbleCon_final,2)/2)
plin_bubble(1:i,:,:,:) = BubbleCon_final(1:i)'./(dRadius*4*pi);
plin_bubble(i+1:end,:,:,:) = BubbleCon_final(i+1:end)'./(dRadius*4*pi);


plin1 = plin;
pdiff  = plin1-plin0;
BubblePos = [find(BubbleLoc_final(1)+TXYstart(2)*dxpar<xpar,1),find((slicesavepos - ypar)<1e-10,1),find(BubbleLoc_final(3)+output.start(4)*dzpar<zpar,1)];

figure(3)
ax = gca;
hold on;
plot(tpar,pdiff(:,BubblePos(1),BubblePos(3)),'--','color','k')
plot(tpar,plin_bubble(:,BubblePos(2),BubblePos(1),BubblePos(3))','-','color',[17 17 17]/255)
ax.YLabel.String = 'Pressure [Pa]'
set(ax,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [?sec]';

% plin_interpx = interp1([xnorm(BubblePos(1)) xnorm(BubblePos(1)-1)],[plin_bubble(:,BubblePos(2),BubblePos(1),BubblePos(3)) plin_bubble(:,BubblePos(2),BubblePos(1)-1,BubblePos(3))]',BubbleLoc_finalN(1));
% plin_interpz = interp1([znorm(BubblePos(3)-1) znorm(BubblePos(3))],[plin_bubble(:,BubblePos(2),BubblePos(1),BubblePos(3)-1) plin_bubble(:,BubblePos(2),BubblePos(1),BubblePos(3))]',BubbleLoc_finalN(3));% plot(BubbleCon_final)
% plot((plin_interpx+plin_interpz)/2)
legend('PDiff','Bubble')
sqrt(sum(sum(sum((pdiff - squeeze(plin_bubble(:,BubblePos(2),:,:))).^2))))