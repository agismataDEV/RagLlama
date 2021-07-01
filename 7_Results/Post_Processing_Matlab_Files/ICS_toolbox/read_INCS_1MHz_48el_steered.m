freq0 = 0.8E6; %fundamental frequency
rootname = 'output_SHI_0MHz8_88el';
dirname = '0MHz8_88el';
beamiterations = [7];
slicesavedim = ['y' 'z'];
slicesavepos = [0 60];
strides = [1 1 1 1];
offsets = [0 0 0 0];

focalplanename = [ slicesavedim(1) int2string_ICS(slicesavepos(1)) ];
transversalplanename = [ slicesavedim(2) int2string_ICS(slicesavepos(2)) ];

% load first beam of focal plane
filename = [ dirname '/' rootname ...
    int2string_ICS(beamiterations(1)-1)...
    '_' focalplanename ...
    int2string_ICS(0) ];

beamdata = load_ICS_slice(filename,strides,offsets);

ppar=squeeze(beamdata.data);
tdimpar=size(ppar,1);
xdimpar=size(ppar,2);
zdimpar=size(ppar,3);
TXYoff=beamdata.TXYoff;
TXYstart=beamdata.start(1:3);
dtpar=beamdata.stepsize(1); % step size dt in s
dxpar=beamdata.stepsize(2) * 1e3; % step size dx in mm
dzpar=beamdata.stepsize(4) * 1e3; % step size dz in mm
zstart= beamdata.start(4)*dxpar; 

clear beamdata

% load data for all higher beams of focal plane
% update the created variables with them
for beam = 1:numel(beamiterations)-1

    filename = [ dirname '/' rootname ...
        int2string_ICS(beamiterations(1+beam)-1)...
        focalplanename ...
        int2string_ICS(beam) ];

    beamdata = load_ICS_slice(filename,strides,offsets);

    ppar(:,:,zdimpar+(1:size(beamdata.data,4)))=squeeze(beamdata.data);
    TXYoff(zdimpar+(1:size(beamdata.data,4)),:)=beamdata.TXYoff+repmat(beamdata.start(1:3)-TXYstart,size(beamdata.data,4),1);
    zdimpar = zdimpar + size(beamdata.data,4);

    clear beamdata
end

% generate axes for focal plane
TXYstarts=TXYstart+min(TXYoff);
tdimpar_glob=tdimpar+max(TXYoff(:,1))-min(TXYoff(:,1));
xdimpar_glob=xdimpar+max(TXYoff(:,2))-min(TXYoff(:,2));
tpar=(0:tdimpar-1)*dtpar;
tpar_glob=(TXYstart(1)+min(TXYoff(:,1))+(0:tdimpar_glob-1))*dtpar;
xpar=(TXYstart(2)+TXYoff(1,2)+(0:xdimpar-1))*dxpar;
xpar_glob=(TXYstart(2)+min(TXYoff(:,2))+(0:xdimpar_glob-1))*dxpar;
zpar = (zstart + (0:zdimpar-1)*dzpar);

% load data of transversal plane
% first find the beam in which the transversal plane was for the filename
[transversalindex,temp] = min(abs((1:zdimpar)*dzpar - slicesavepos(2)));
transversalbeam = floor(transversalindex/ceil(zdimpar/numel(beamiterations)));
filename = [ dirname '/' rootname ...
    int2string_ICS(beamiterations(1+transversalbeam)-1)...
    '_' transversalplanename ];

beamdata = load_ICS_slice(filename,[1 1 1 1],[0 0 0 0]);

size_trans=size(beamdata.data);
ppar_trans=zeros(size_trans(1),size_trans(2),2*size_trans(3)-1);
ppar_trans(:,:,size_trans(3):2*size_trans(3)-1)=squeeze(beamdata.data);
ppar_trans(:,:,size_trans(3)-1:-1:1)=ppar_trans(:,:,size_trans(3)+1:2*size_trans(3)-1);
tdimpar_tr=size(ppar_trans,1);
xdimpar_tr=size(ppar_trans,2);
ydimpar_tr=size(ppar_trans,3);
TXYstart_tr=beamdata.start(1:3);
TXYstart_tr(3)=-size_trans(3)+1;
dtpar_tr=beamdata.stepsize(1);
dxpar_tr=beamdata.stepsize(2) * 1e3;
dypar_tr=beamdata.stepsize(2) * 1e3;
z_tr=beamdata.start(4)*beamdata.stepsize(4);

clear beamdata

% generate axes for transversal plane
xpar_trans=(TXYstart_tr(2)+(0:xdimpar_tr-1))*dxpar_tr;
ypar=(TXYstart_tr(3)+(0:ydimpar_tr-1))*dxpar_tr;

% echo the MI - approx (coarse discretization)
mechanical_index = max(max(max(-ppar)))/1e6/sqrt(freq0/1e6);
disp(sprintf('Mechanical Index: %0.4g',mechanical_index))

% store the pressure at the focus, (x,z) = (0,60) mm
eval(['p_' dirname ' = ppar_trans(:,floor(xdimpar_tr/2)+1,floor(ydimpar_tr/2)+1);']); 
eval(['save ''p_' dirname ''' p_' dirname]);

% generate fundamental, higher harmonic and superharmonic profiles
pparfund=nan(xdimpar_glob,zdimpar);
ppar2nd=nan(xdimpar_glob,zdimpar);
% ppar3rd=nan(xdimpar_glob,zdimpar);
% ppar4th=nan(xdimpar_glob,zdimpar);
% ppar5th=nan(xdimpar_glob,zdimpar);
pparSH=nan(xdimpar_glob,zdimpar);
for i=1:zdimpar
    pparfund((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(ppar(:,:,i)),dtpar,.7*freq0,1.3*freq0));
    ppar2nd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(ppar(:,:,i)),dtpar,1.7*freq0,2.3*freq0));
%     ppar3rd((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(ppar(:,:,i)),dtpar,2.7*freq0,3.3*freq0));
%     ppar4th((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(ppar(:,:,i)),dtpar,3.7*freq0,4.3*freq0));
%     ppar5th((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(ppar(:,:,i)),dtpar,4.7*freq0,5.3*freq0));
    pparSH((1:xdimpar)+TXYoff(i,2)-min(TXYoff(:,2)),i)=20*log10(maxhilbert(squeeze(ppar(:,:,i)),dtpar,2.7*freq0,5.3*freq0));
end
maxfund=max(max(pparfund));
% 
% save 'SHI_10MHz_1o2_newtx_0deg_profiles' pparfund ppar2nd ppar3rd ppar4th ppar5th pparSH maxfund ppar zpar xpar_glob tpar
% 
% ppartrfund=nan(xdimpar_tr,ydimpar_tr);
% ppartr2nd=nan(xdimpar_tr,ydimpar_tr);
% % ppartr3rd=nan(xdimpar_tr,ydimpar_tr);
% % ppartr4th=nan(xdimpar_tr,ydimpar_tr);
% % ppartr5th=nan(xdimpar_tr,ydimpar_tr);
% ppartrSH=nan(xdimpar_tr,ydimpar_tr);
% for i=1:ydimpar
%     ppartrfund(:,i)=20*log10(maxhilbert(squeeze(ppar_trans(:,:,i)),dtpar_tr,.7*freq0,1.3*freq0));
%     ppartr2nd(:,i)=20*log10(maxhilbert(squeeze(ppar_trans(:,:,i)),dtpar_tr,1.7*freq0,2.3*freq0));
% %     ppartr3rd(:,i)=20*log10(maxhilbert(squeeze(ppar_trans(:,:,i)),dtpar_tr,2.7*freq0,3.3*freq0));
% %     ppartr4th(:,i)=20*log10(maxhilbert(squeeze(ppar_trans(:,:,i)),dtpar_tr,3.7*freq0,4.3*freq0));
% %     ppartr5th(:,i)=20*log10(maxhilbert(squeeze(ppar_trans(:,:,i)),dtpar_tr,4.7*freq0,5.3*freq0));
%     ppartrSH(:,i)=20*log10(maxhilbert(squeeze(ppar_trans(:,:,i)),dtpar_tr,2.7*freq0,5.3*freq0));
% end
% maxfund_tr=max(max(ppartrfund));
% % 
% % save 'SHI_10MHz_1o2_newtx_0deg_tr_profiles' ppartrfund ppartr2nd ppartr3rd ppartr4th ppartr5th ppartrSH maxfund_tr ppar_trans z_trans xpar_trans ypar tpar
% % 
return

figure(1),surf(zpar,xpars,pparsXfund-maxfund);shading flat;
axis equal;axis tight
axis([0 150 -20 20 -400 0])
view(0,90);
caxis([-40 0])
% colorbar('YTick',[-40 -30 -20 -10 0],...
%     'YTickLabel',{'-40dB','-30dB','-20dB','-10dB','   0dB'});
% xlabel('z (mm)')
% ylabel('x (mm)')
% title('Main Beam, F0')
set(gca,'XTick',-1)
set(gca,'YTick',-21)
set(gca,'Layer','top')
set(gca,'Box','on')

bewaarsurf(1,'0deg_profile_F0',16,150)
% bewaarsurf(1,'0deg_profile_F0_LPR',16,300)

figure(2),surf(zpar,xpars,pparsX2nd-maxfund);shading flat;
axis equal;axis tight
axis([0 150 -20 20 -400 0])
caxis(max(max(pparsX2nd-maxfund))+[-40 0]);
view(0,90);
% colorbar('YTick',[-60 -50 -40 -30 -20],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB','-20dB'});
% xlabel('z (mm)')
% ylabel('x (mm)')
% title('Main Beam, 2H')
set(gca,'XTick',-1)
set(gca,'YTick',-21)
set(gca,'Layer','top')
set(gca,'Box','on')

% bewaarsurf(2,'0deg_profile_TH',16,150)
bewaarsurf(2,'0deg_profile_TH_LPR',16,300)

figure(3),surf(zpar,xpars,pparsXSH-maxfund);shading flat;
axis equal;axis tight
axis([0 150 -20 20 -400 0])
view(0,90);
caxis(max(max(pparsXSH-maxfund))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
% xlabel('z (mm)')
% ylabel('x (mm)')
% title('Main Beam, SH')
set(gca,'XTick',-1)
set(gca,'YTick',-21)
set(gca,'Layer','top')
set(gca,'Box','on')

% bewaarsurf(3,'0deg_profile_SH',16,150)
bewaarsurf(3,'0deg_profile_SH_LPR',16,300)

figure(7),surf(ypar,xpar_trans,pparstrfund-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis([-40 0]);
colorbar('YTick',[-40 -30 -20 -10 0],...
    'YTickLabel',{'-40dB','-30dB','-20dB','-10dB','   0dB'});
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=42.5mm, F0')

bewaarsurf(7,'0deg_transversal_F0',16,150)

figure(8),surf(ypar,xpar_trans,pparstr2nd-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparstr2nd-maxfund))+[-40 0]);
colorbar('YTick',[-50 -40 -30 -20],...
    'YTickLabel',{'-50dB','-40dB','-30dB','-20dB'});
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=42.5mm, 2H')

bewaarsurf(8,'0deg_transversal_2H',16,150)

figure(9),surf(ypar,xpar_trans,pparstrSH-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparstrSH-maxfund))+[-40 0]);
colorbar('YTick',[-60 -50 -40 -30 -20],...
    'YTickLabel',{'-60dB','-50dB','-40dB','-30dB','-20dB'});
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=42.5mm, SH')

bewaarsurf(9,'0deg_transversal_SH',16,150)

return

% load 'SHI_10MHz_1o2_newtx_0deg_profiles'
% load 'SHI_10MHz_1o2_newtx_0deg_tr_profiles'

figure(1),surf(zpar,xpars,pparsXfund-maxfund);shading flat;view(0,90);
axis equal;axis tight
caxis([-40 0]);
colorbar('YTick',[-40 -30 -20 -10 0],...
    'YTickLabel',{'-40dB','-30dB','-20dB','-10dB','   0dB'});
xlabel('z (mm)')
ylabel('x (mm)')
title('Main Beam, F0')

bewaarsurf(1,'ppar10MHz_1o2_newtx_0deg_F0',16,300)

figure(2),surf(zpar,xpars,pparsX2nd-maxfund);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparsX2nd-maxfund))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('z (mm)')
ylabel('x (mm)')
title('Main Beam, H2')

bewaarsurf(2,'ppar10MHz_1o2_newtx_0deg_2H',16,300)

figure(3),surf(zpar,xpars,pparsX3rd-maxfund);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparsX3rd-maxfund))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('z (mm)')
ylabel('x (mm)')
title('Main Beam, H3')

bewaarsurf(3,'ppar10MHz_1o2_newtx_0deg_3H',16,300)

figure(4),surf(zpar,xpars,pparsX4th-maxfund);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparsX4th-maxfund))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('z (mm)')
ylabel('x (mm)')
title('Main Beam, H4')

bewaarsurf(4,'ppar10MHz_1o2_newtx_0deg_4H',16,300)

figure(5),surf(zpar,xpars,pparsX5th-maxfund);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparsX5th-maxfund))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('z (mm)')
ylabel('x (mm)')
title('Main Beam, H5')

bewaarsurf(5,'ppar10MHz_1o2_newtx_0deg_5H',16,300)

figure(6),plot(tpar*1e6,squeeze(ppar(:,71,221))/1e6)
axis([3 12 -1.7 1.7])
xlabel('t (mus)')
ylabel('p (MPa)')
title('Pulse at the focus (x,z)=(0,60) mm')

bewaar(6,'ppar10MHz_1o2_newtx_0deg_pulse',14)

figure(7),surf(ypar,xpar_trans,pparstrfund-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis([-40 0]);
colorbar('YTick',[-40 -30 -20 -10 0],...
    'YTickLabel',{'-40dB','-30dB','-20dB','-10dB','   0dB'});
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=60mm, F0')

bewaarsurf(7,'ppar10MHz_1o2_newtx_0deg_tr_F0',16,300)

figure(8),surf(ypar,xpar_trans,pparstr2nd-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparstr2nd-maxfund_tr))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=60mm, H2')

bewaarsurf(8,'ppar10MHz_1o2_newtx_0deg_tr_2H',16,300)

figure(9),surf(ypar,xpar_trans,pparstr3rd-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparstr3rd-maxfund_tr))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=60mm, H3')

bewaarsurf(9,'ppar10MHz_1o2_newtx_0deg_tr_3H',16,300)

figure(10),surf(ypar,xpar_trans,pparstr4th-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparstr4th-maxfund_tr))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=60mm, H4')

bewaarsurf(10,'ppar10MHz_1o2_newtx_0deg_tr_4H',16,300)

figure(11),surf(ypar,xpar_trans,pparstr5th-maxfund_tr);shading flat;view(0,90);
axis equal;axis tight
caxis(max(max(pparstr5th-maxfund_tr))+[-40 0]);
% colorbar('YTick',[-60 -50 -40 -30],...
%     'YTickLabel',{'-60dB','-50dB','-40dB','-30dB'});
colorbar
xlabel('y (mm)')
ylabel('x (mm)')
title('Main Beam, z=60mm, H5')

bewaarsurf(11,'ppar10MHz_1o2_newtx_0deg_tr_5H',16,300)
