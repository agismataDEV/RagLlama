
load('../Results_for_separate_Variables/phi.mat')
%%
clear an
an.medium.freq0         = 15E6;
an.medium.Fnyq          = domain.Fnyq;
an.medium.rho0          = 1060;
an.medium.c0            = 1480;
an.medium.kappa0        = 2*pi*an.medium.freq0/an.medium.c0;

an.domain.dtpar         = 1/(2*an.medium.Fnyq*an.medium.freq0);
an.domain.dims          = size(plin.data);
an.domain.tdimpar       = an.domain.dims(1)
an.domain.tpar          = [0:an.domain.dims(1)-1]*an.domain.dtpar;

an.domain.dxpar         = an.domain.dtpar*an.medium.c0;
an.domain.zpar          = ([1:an.domain.dims(3)] - domain.start(3))*an.domain.dxpar;
an.domain.xstart        =  -an.domain.dims(2)/2;
an.domain.xpar          = (an.domain.xstart+(0:an.domain.dims(2)-1))*an.domain.dxpar;%an.domain.xpar = domain.par{1};

an.domain.beam          = 0;
an.file.dirname         = file.dirname;
% an.file.dirname         = '../David_210deg_x_Lagrangian_Only';
an.file.focalplanename  = [ 'y' int2string_ICS(2) ];

Fs          = 1/an.domain.dtpar;  dF = Fs/an.domain.tdimpar;
if mod(an.domain.tdimpar,2) == 0 % Nt is even
    f = -Fs/2 : dF : Fs/2-dF;
else            % Nt is odd
    f = [sort(-1*(dF:dF:Fs/2)) (0:dF:Fs/2)];
end
an.omega       = 2*pi*f;
an.omega       = [an.omega(an.omega>=0) an.omega(an.omega<0)];
an.theta       = 21;
%% Symbolics
an.depth                = 258;
an.lat                  = 169;

syms t z x
an.blackman_L             = 6;
an.blackman_delay         = 6;
an.blwin                  = 0.35869 - 0.48829 * cos(2*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) + 0.14128*cos(4*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) - 0.01168*sin(6*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L);
an.sign                   = (sign(t-an.blackman_delay/an.medium.freq0)-sign(t-an.blackman_delay/an.medium.freq0-an.blackman_L/an.medium.freq0))/2;

an.p_expr                 = 4E4*sin(2*pi*an.medium.freq0*t-6/an.medium.freq0 - an.medium.kappa0 * z)  * an.blwin * an.sign;

an.gauss_dl                  = 6.0/an.medium.freq0;
an.gauss_win                 = 3.0/an.medium.freq0;
w                            = 2 * pi * an.medium.freq0;

A                         = ((t-an.gauss_dl)/(an.gauss_win/2));
an.p_expr                 = 1E5*exp( - A.^2).* sin(w*(t-gauss_dl) - an.medium.kappa0*sind(180 - theta)*x - an.medium.kappa0*cosd(180-theta)*z);
% an.p_expr                 = 1E5*exp( - A.^2).*(sin(w*(t-gauss_dl) - an.medium.kappa0*sind(theta)*x + an.medium.kappa0*cosd(theta)*z)*-A/(gauss_win/2*an.medium.freq0) + cos(w*(t-gauss_dl)- an.medium.kappa0*sind(theta)*x + an.medium.kappa0*cosd(theta)*z));
an.p_expr                 = 1E5*(sin(w*(t-gauss_dl) - an.medium.kappa0*sind(180 - theta)*x - an.medium.kappa0*cosd(180 - theta)*z) + sin(w*(t-gauss_dl) + an.medium.kappa0*sind(180-theta)*x- an.medium.kappa0*cosd(180-theta)*z));
% an.p_expr                 = 1E5*(sin(w*(t-gauss_dl) - an.medium.kappa0*sind(180 - theta)*x - an.medium.kappa0*cosd(180 - theta)*z));
an.phi_expr               = int(-an.p_expr,t)/an.medium.rho0*an.medium.freq0;

%UN-SIMPLIFIED EXPRESSION VARIABLES
an.nabla_z_sq_expr        = (diff(an.phi_expr,z))^2/an.medium.freq0^2;
an.nabla_x_sq_expr        = (diff(an.phi_expr,x))^2/an.medium.freq0^2;
an.nabla_t_sq_expr        = (an.p_expr/an.medium.c0/an.medium.rho0)^2;%(diff(an.phi_expr,t)/an.medium.c0/an.medium.freq0)^2;

%SIMPLIFIED EXPRESSION VARIABLES
an.phi_sq_expr            = an.phi_expr.^2;
an.lapl_t_expr            = diff(diff(an.phi_sq_expr,t),t)/an.medium.c0^2;
an.lapl_x_expr            = diff(diff(an.phi_sq_expr,z),z);

%%
% GET THE VALUES OF EACH VARIABLE FOR EACH TIMEPOINT
an.p                    = double(subs(subs(an.p_expr           , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.phi                  = double(subs(subs(an.phi_expr         , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.phi_sq               = double(subs(subs(an.phi_sq_expr      , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.lapl_t               = double(subs(subs(an.lapl_t_expr      , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.lapl_x               = double(subs(subs(an.lapl_x_expr      , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.nabla_x_sq           = double(subs(subs(an.nabla_x_sq_expr  , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.nabla_z_sq           = double(subs(subs(an.nabla_z_sq_expr  , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.nabla_t_sq           = double(subs(subs(an.nabla_t_sq_expr  , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
%%
[an.nabla_x_sq,an] = Compare_Signals(an,plin,pnl,domain,dslice,file,an.nabla_x_sq);
[an.nabla_z_sq,an] = Compare_Signals(an,plin,pnl,domain,dslice,file,an.nabla_z_sq);
[an.nabla_t_sq,an] = Compare_Signals(an,plin,pnl,domain,dslice,file,an.nabla_t_sq);
close all
%%
[an.p,an] = Compare_Signals(an,plin,pnl,domain,dslice,file,an.p);
%%  ====================================================== P ======================================================

pnl.an.filename     = [an.file.dirname '/' 'TESTNeumann' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.p]              = FieldLoad(an.filename);
Frequency(domain, sim.p, an.p , an.depth, sim.title);

%%  ====================================================== PHI ======================================================

pnl.an.filename     = [an.file.dirname '/' 'phi_fourier' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.phi]            = FieldLoad(an.filename);
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(phi, sim.phi, an.phi, an.depth, '\phi as a function of time' );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMPLIFIED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ====================================================== PHI_SQ ======================================================

an.filename        = [an.file.dirname '/' 'phi_sq_fourier' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.phi_sq]        = FieldLoad(an.filename);
Frequency(domain, sim.phi_sq,real(an.phi_sq),an.depth, sim.title);

%% ====================================================== Time Derivative of (phi)^2 ============================================================

an.filename        = [an.file.dirname '/' 'Time_Deriv' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.lapl_t]        = FieldLoad(an.filename);
% phi.td_an = real(double(subs(an.lapl_t_expr,t,an.domain.tpar)));
% phi.td_an(isnan(phi.td_an))=0;
Frequency(domain, sim.lapl_t,real(an.lapl_t)/an.medium.freq0^2,an.depth, sim.title);
% '$\frac{1}{c_0^2}\frac{\vartheta^2 \phi^2}{\vartheta t^2} [m^2/sec^2]$','Interpreter','latex','FontSize',40)

%%  ====================================================== Nabla^2 of (phi)^2 ======================================================

an.filename        = [an.file.dirname '/' 'Z_Component' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.lapl_x]        = FieldLoad(an.filename);
Frequency(domain,sim.lapl_x,real(an.lapl_x),an.depth, sim.title);

%% ====================================================== Lagrangian ======================================================

an.filename        = [an.file.dirname '/' 'Lagrangian' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.lagrangian]    = FieldLoad(an.filename);
an.lagrangian = real(an.lapl_t/an.medium.freq0^2 - an.lapl_x);
an.lagrangian(isnan(an.lagrangian))=0;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(domain, sim.lagrangian,an.lagrangian,an.depth, sim.title);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ======================================================(Time Derivative of phi)^2 ======================================================

an.filename        = [an.file.dirname '/' 'Time_Deriv' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.nabla_t_sq]    = FieldLoad(an.filename);
% an.nabla_t_sq = real(double(subs(an.nabla_t_sq_expr,t,an.domain.tpar)));
% an.nabla_t_sq(isnan(an.nabla_t_sq))=0;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(domain, sim.nabla_t_sq,an.nabla_t_sq,an.depth, sim.title);

%% ====================================================== (Nabla of phi ) ^2 ======================================================

an.filename        = [an.file.dirname '/' 'Z_Component' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.nabla_x_sq]    = FieldLoad(an.filename);
Frequency(domain, sim.nabla_x_sq,real(an.nabla_x_sq),an.depth, sim.title);

%% ====================================================== Lagrangian ======================================================

an.filename        = [an.file.dirname '/' 'Lagrangian' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.lagrangian]    = FieldLoad(an.filename);
an.lagrangian       = real(an.nabla_t_sq - an.nabla_x_sq)*an.medium.rho0/2;
an.lagrangian(isnan(an.lagrangian))=0;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.lagrangian, an.lagrangian', an.depth, sim.title);

%% PLOT RANDOM variables
sim.title = ['Pressure due to Lagrangian @ 0 deg, Focused'];

an.filename        = [an.file.dirname '/' 'LagrangianPressure' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.rnd_var]       = FieldLoad(an.filename);
% an.rnd_var          = real(an.nabla_t_sq - an.nabla_x_sq)*an.medium.rho0/2;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.rnd_var, an.rnd_var', an.depth, sim.title);
GENERATE_MAP(an, sim.rnd_var, sim.title);

%% ------------------------------------------------------ GENERAL FUNCTIONS -------------------------------------------------------------------------
function [fourier] = FieldLoad(filename)

output_nonl         = load_ICS_slice(filename);
fourier             = squeeze(output_nonl.data);

end

function  [f_p,S_sim,S_an] = Frequency(an, sim_var, an_var ,depth, title_var)

F_s = 1/(2*an.domain.dtpar);
% Calculate the frequency spectrum based on a sampling frequency and the pressure values
Nfft_p = 2^nextpow2(2*length(an.domain.tpar));                      % Zero-padding to improve fft performance until length equal to 2N
Y_phi_fourier = abs(fft(sim_var(:,an.lat,depth),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
f_p = linspace(0,1,Nfft_p/2+1)*F_s;                                                           % [Hz] Frequency sample Number
% Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
S_sim = 20*log10(Y_phi_fourier(Np));

Y_phi_an = abs(fft(an_var,Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
S_an = 20*log10(Y_phi_an(Np));

figure('WindowState','maximized'); hold on; grid on;
title('Frequency spectrum , grad phi','Interpreter','latex')
plot(f_p*1E-6,S_sim-max(S_sim),'--','Color','blue','LineWidth',2)
plot(f_p*1E-6,S_an-max(S_an),'k-','LineWidth',2)
legend('INCS','Analytical','Interpreter','latex')
xlabel('Frequency [MHz]','Interpreter','latex')
ylabel('Amplitude [dB]','Interpreter','latex')
ylim([-80 0])
title([title_var , ' @ (X,Y,Z) =(', num2str(an.domain.xpar(an.lat)*1E3),', ', num2str(0),', ', num2str(an.domain.zpar(an.depth)*1E3),') [mm]'],'Interpreter','latex')

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',15)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)

an.OSFactor = 5;
an.domain.tpar_OS = (0:size(an_var,1)*an.OSFactor-1)*an.domain.dtpar/an.OSFactor;

figure('WindowState','maximized');hold on; grid on;
title('grad phi computed via different methods','Interpreter','latex')
plot(an.domain.tpar_OS *1E6, resampleSINC(sim_var(:,an.lat,depth)',an.OSFactor),'--','Color','blue','LineWidth',2)
plot(an.domain.tpar_OS *1E6, resampleSINC(an_var',an.OSFactor),'k-','LineWidth',2)

legend('INCS','Analytical','Interpreter','latex')
xlabel('Time [$\mu$s]','Interpreter','latex')
ylabel('phi [m^2/sec]','Interpreter','latex')
title([title_var , ' @ (X,Y,Z) =(', num2str(an.domain.xpar(an.lat)*1E3),', ', num2str(0),', ', num2str(an.domain.zpar(an.depth)*1E3),') [mm]'],'Interpreter','latex')

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',15)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)
end

function GENERATE_MAP(an, sim_var, title_var)

%==================================================== UNFOCUSED
figure('WindowState','maximized');
length_z = 1;
imagesc(an.domain.zpar(length_z:an.domain.dims(3) -length_z)*1E3,an.domain.xpar*1E3,squeeze(max(abs(sim_var(:,:,length_z:end-length_z)))));

xlabel('Z [mm]','Interpreter','latex')
ylabel('X [mm]','Interpreter','latex' )
title(title_var,'Interpreter','latex')
c=colorbar;
c.Label.String = 'Lagrangian Density [kg/(m$\cdot$s)]';
c.Label.Interpreter='latex';
c.TickLabelInterpreter='latex';

ax =gca;
ax.TickLabelInterpreter = 'latex';
file.plot_colour = 'viridis';xAM_Colormaps(file)
set(gca,'FontSize',20)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)
%==================================================== FOCUSED
figure('WindowState','maximized')
start_z = 100;
end_z = an.domain.dims(3)-40;
imagesc(an.domain.zpar(start_z:end_z)*1E3,an.domain.xpar*1E3,squeeze(max(abs(sim_var(:,:,start_z:end_z)))));

xlabel('Z [mm]','Interpreter','latex' )
ylabel('X [mm]','Interpreter','latex' )
title([title_var, ', Zoom-In'],'Interpreter','latex')

c=colorbar;
c.Label.String = 'Lagrangian Density [kg/(m$\cdot$s)]';
c.Label.Interpreter='latex';
c.TickLabelInterpreter='latex';
file.plot_colour = 'viridis';xAM_Colormaps(file)

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)
end

function plot_time_frame(an,pnl,domain,dslice,file,frame)
%%
i_end = an.domain.tdimpar;
plot_value = (20*log10(abs(pnl.xwave_data)));
dBVALUES = 0.3*max(max(max(plot_value)));

%% ========================= Initialize first frame=========================
% This is done to speedup the process of making images , the most time is lost in making the plots
num_frames = i_end;
F(num_frames) = struct('cdata',[],'colormap',[]);

f2 = figure('WindowState','maximized');
ax = gca;
p=imagesc(ax,domain.par{domain.dimval(3)},domain.par{domain.dimval(2)},squeeze(plot_value(1,:,:)) );
hold on;
caxis(ax,[max(max(max(plot_value)))-dBVALUES max(max(max(plot_value)))])
title_txt = 'Left-Aperture Pressure Field, $t_i$ = ';
title(ax,[title_txt,num2str(1)],'Interpreter','latex')
xlabel(ax,dslice.xlabel,'Interpreter','latex')
ylabel(ax,dslice.ylabel,'Interpreter','latex')

xlim(ax,[min(domain.par{domain.dimval(3)}) max(domain.par{domain.dimval(3)})])
ylim(ax,[min(domain.par{domain.dimval(2)}) max(domain.par{domain.dimval(2)})])

xAM_Colormaps(file)
c = colorbar;
c.Label.String = 'Pressure [dB]';
c.Label.Interpreter ='latex';
c.TickLabelInterpreter = 'latex';

set(ax,'YDir','normal')
if dslice.savedim(dslice.num) =='z'; set(ax,'XDir','reverse');set(ax,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
axis image;
set(gcf, 'Color', 'w');
set(gca,'FontSize',18)
ax.TickLabelInterpreter = 'latex';
set(findall(gca,'-property','FontSize'),'FontSize',30)

for i=frame:frame % Campure a specific time frame and plot it in a figure
    set(ax.Title,'String',[title_txt,num2str(i)],'Interpreter','latex')
    set(p,'XData',domain.par{domain.dimval(3)},'YData',domain.par{domain.dimval(2)},'CData',squeeze(plot_value(i,:,:)))
    drawnow
end
plot([domain.par{3}(1) domain.par{3}(end)],[0 0],'--','Color','white','LineWidth',2)
%Capture the indexes of the max values in this plot
[max_val,max_idx] = max(squeeze(plot_value(i,:,:))');

% Plot the values in a figure
f3 = figure('WindowState','maximized');
plot(domain.par{3}(max_idx),domain.par{1},'LineWidth',2)

xlabel(dslice.xlabel,'Interpreter','latex')
ylabel(dslice.ylabel,'Interpreter','latex')
axis image;
xlim([domain.par{3}(1) domain.par{3}(end)])

% Add 2 points in the same figure
hold on;
P_ZX = [domain.par{3}(max_idx) ; domain.par{1}]';
P1 = P_ZX(160,:);
P2 = P_ZX(191,:);
plot(P1(:,1),P1(:,2),'ro','LineWidth',2','MarkerSize',5)
plot(P2(:,1),P2(:,2),'ro','LineWidth',2','MarkerSize',5)

%Find many differences
P1_many = P_ZX(160,:);
P2_many = P_ZX(191,:);
Angle_many=90-atand((P2_many(:,2)-P1_many(:,2))./(P2_many(:,1)-P1_many(:,1)));
title(['Max pressure location in X-Z plane, ' num2str(domain.angle) ' deg, @ Y = 0 mm '],'Interpreter','latex')
%     title(['Average angle between the 2 points, ',num2str(mean(Angle_many)) ,' deg'])
legend(['Angle between 2 points: ',num2str(mean(Angle_many)) ,' deg'],'Interpreter','latex')
plot([domain.par{3}(1) domain.par{3}(end)],[0 0],'--','Color','black','LineWidth',2)

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(findall(gca,'-property','FontSize'),'FontSize',40)

set(gcf, 'Color', 'w');
set(gca,'FontSize',18)
end

function [an_press,an] = Compare_Signals(an,plin,pnl,domain,dslice,file,var)
%%
[x_el,arr,td_x,apod_x,apod_right,apod_left] = xAM_Generate_TD_APOD(an.theta);
dist_cross_x_axis = an.domain.zpar(an.depth);
if (an.theta>0)
    element_n = 13;
    element_x = abs(x_el(element_n))*1E-3; % m
    dist_cross_x_axis = element_x*tand(90-an.theta); % m
    [val,an.depth]= min(abs(an.domain.zpar-dist_cross_x_axis));
else
    td_x = zeros(32,1);element_n = 20;element_x =0;
end
travel_dist_xwave = dist_cross_x_axis*cosd(an.theta) + element_x*sind(an.theta) + domain.start(3)*an.domain.dxpar - an.domain.dxpar; % m
an_td = travel_dist_xwave/an.medium.c0 + td_x(32 - element_n +1)*1E-6 + an.domain.dtpar; % There is some small difference because the dist that cross x is not always a gridpoint in the computational domain
frame = round((an_td + an.gauss_dl)/an.domain.dtpar );
plot_time_frame(an,pnl,domain,dslice,file,frame)

window = an.gauss_dl * an.medium.freq0 * 2 * an.medium.Fnyq + an.depth/an.medium.c0;
len_window  = (an.gauss_win + an.gauss_dl)  * an.medium.freq0 * 2 * an.medium.Fnyq;

rect_win  = [zeros(round(an_td/an.domain.dtpar),1); rectwin(len_window); zeros(round(an.domain.tdimpar-an_td/an.domain.dtpar - len_window),1)];
an_press  = real(ifft(fft(var).*exp(-1i*an.omega*an_td))).* rect_win';

Frequency(an, plin.data, an_press' ,an.depth, ['Signal Comparison']);

figure('WindowState','maximized');
hold on;plot(an.domain.xpar*1E3,max(abs(plin.data(:,:,an.depth)))*1E-3,'LineWidth',2)

ax = gca;
xlabel('X [mm]','Interpreter','latex')
ylabel('Pressure [kPa]','Interpreter','latex')
title(['Lateral profile @ Z =', num2str(an.domain.zpar(an.depth)*1E3),' [mm]'],'Interpreter','latex')

ax.TickLabelInterpreter = 'latex';
set(findall(gca,'-property','FontSize'),'FontSize',40)
set(gcf,'Color','white')
end