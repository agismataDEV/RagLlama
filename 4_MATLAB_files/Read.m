%% L.D. 10-05-2012 - This script reads and displays the results contained in the HDF5 data
clear all
% close all
clc

%% INITIALIZATION
IterationNumber     = 5;                                                   % iteration number considered
name                = '100_35FNYQ_NFBU_FA/TESTNeumann'                                        % file name
f0                  = 1.0;                                                  % fundamental frequency
f0                  = f0*1000000; % freq in Hz
ii                  = 0;

%% LINEAR DATA
first       = [name '_000_y_001_00' num2str(ii)];
    
output      = load_parnac_slice(first);
plin        = squeeze(output.data);
tdimpar     = size(plin,1);
xdimpar     = size(plin,2);
zdimpar     = size(plin,3);
ydimpar     = size(plin,3);
TXYoff      = output.TXYoff;
TXYstart    = output.start(1:3);
dtpar       = output.stepsize(1);
dxpar       = output.stepsize(2);
dypar       = output.stepsize(3);
dzpar       = output.stepsize(4);
zstart      = output.start(4)*dzpar;
xsource     = (output.start(2)+(0:size(plin,2)-1))*dxpar*1000;
ysource     = (output.start(3)+(0:size(plin,3)-1))*dypar*1000;
  
%% NONLINEAR DATA
third       = [name '_00' num2str(IterationNumber) '_y_001_00' num2str(ii)];

if (IterationNumber>9)&&(IterationNumber<100)
    third       = [name '_0' num2str(IterationNumber) '_y_001_00' num2str(ii)];
else if IterationNumber>99
        third       = [name '_' num2str(IterationNumber) '_y_001_00' num2str(ii)];
    end
end
 
output      = load_parnac_slice(third); 
pnonl       = squeeze(output.data); 
xdimpars    = xdimpar+max(TXYoff(:,2))-min(TXYoff(:,2)); 
ydimpars    = ydimpar+max(TXYoff(:,3))-min(TXYoff(:,3));
ypar        = (TXYstart(3)+TXYoff(1,3)+(0:ydimpar-1))*dypar*1000; 
ypars       = (TXYstart(3)+min(TXYoff(:,3))+(0:ydimpars-1))*dypar*1000;
tpar        = (0:tdimpar-1)*dtpar; 
xpar        = (TXYstart(2)+TXYoff(1,2)+(0:xdimpar-1))*dxpar*1000;
xpars       = (TXYstart(2)+min(TXYoff(:,2))+(0:xdimpars-1))*dxpar*1000;
zpar        = (zstart + (0:zdimpar-1)*dxpar)*1000;

Dimensions=size(pnonl)
sizeP=size(pnonl);
     
%% CREATING MATRICES 
plin_fund = nan(xdimpars,zdimpar);
pnl_fund  = nan(xdimpars,zdimpar);
pnl_2nd   = nan(xdimpars,zdimpar);
pnl_3nd   = nan(xdimpars,zdimpar);
pnl_4nd   = nan(xdimpars,zdimpar);

   
  for i=1:zdimpar
      
    plin_fund((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i) = 20*log10(maxhilbert(squeeze(plin(:,:,i)),dtpar,.6*f0,1.4*f0));
      pnl_fund((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)  = 20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,.6*f0,1.4*f0));   
      pnl_2nd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)   = 20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,1.6*f0,2.4*f0));
      pnl_3nd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)   = 20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,2.6*f0,3.4*f0));
%       pnl_4nd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)   = 20*log10(maxhilbert(squeeze(pnonl(:,:,i)),dtpar,3.6*f0,4.4*f0));
      
  end

%% PLOTS
  dimensions = size(pnl_fund)
  figure(1)
  hold on
  plot(zpar-zpar(1),plin_fund( round(dimensions(1)/2),:),'k')
  plot(zpar-zpar(1),pnl_fund( round(dimensions(1)/2),:),'b')
  plot(zpar-zpar(1),pnl_2nd(  round(dimensions(1)/2),:),'r')
  plot(zpar-zpar(1),pnl_3nd(  round(dimensions(1)/2),:),'y')
  plot(zpar-zpar(1),pnl_4nd(  round(dimensions(1)/2),:),'g')
  grid on
  xlabel('z [m]')
  ylabel('Pressure [dB re 1 Pa] ')  
  legend('Linear Field','F0', '2H', '3H', '4H')
  axis([0 20 50 140 ])
   
  
  dBVALUES=15;
  figure(2)
  subplot(4,1,1)
  imagesc(zpar,xpars,pnl_fund)
  caxis([max(max(pnl_fund))-dBVALUES max(max(pnl_fund))])
  axis image
  title('F0')
  colorbar
  figure(2)
  subplot(4,1,2)
  imagesc(zpar,xpars,pnl_2nd)
  caxis([max(max(pnl_2nd))-dBVALUES max(max(pnl_2nd))])
  title('2H')
  axis image
  colorbar
  figure(2)
  subplot(4,1,3)
  imagesc(zpar,xpars,pnl_3nd)
  caxis([max(max(pnl_3nd))-dBVALUES max(max(pnl_3nd))])
  title('3H')
  axis image
  colorbar
  figure(2)
  subplot(4,1,4)
  imagesc(zpar,xpars,pnl_4nd)
  caxis([max(max(pnl_4nd))-dBVALUES max(max(pnl_4nd))])
  title('3H')
  axis image
  colorbar
  
  %%
figure(3)
ax = gca;
hold on;
plot(tpar,squeeze(pnonl(:,32,1)))
set(ax,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [?sec]';
set(ax,'YTickLabel',ax.YTick*1e-3) ; ax.YLabel.String = 'Pressure [kPa]';
legend('Time Signature')


