
load('../Results_for_separate_Variables/phi.mat')
%%
clear an
an.medium.freq0         = 15E6;
an.medium.Fnyq          = 3;
an.medium.rho0          = 1060;
an.medium.c0            = 1480;
an.medium.kappa0        = 2*pi*an.medium.freq0/an.medium.c0;

an.domain.dtpar         = 1/(2*an.medium.Fnyq*an.medium.freq0);
an.domain.dims          = size(sim.rnd_var);
an.domain.tpar          = [0:an.domain.dims(1)-1]*an.domain.dtpar;

an.domain.dxpar         = an.domain.dtpar*an.medium.c0;
an.domain.zpar          = [0:an.domain.dims(3)-1]*an.domain.dxpar;
an.domain.xstart        =  -an.domain.dims(2)/2;
an.domain.xpar          = (an.domain.xstart+(0:an.domain.dims(2)-1))*an.domain.dxpar;%an.domain.xpar = domain.par{1};
an.depth                = 200;

an.domain.beam          = 0;
an.file.dirname         = file.dirname;
% an.file.dirname         = '../David_210deg_x_Lagrangian_Only';
an.file.focalplanename  = [ 'y' int2string_ICS(2) ];
%% Symbolics
syms t z 
an.blackman_L             = 6;
an.blackman_delay         = 6; 
an.blwin                  = 0.35869 - 0.48829 * cos(2*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) + 0.14128*cos(4*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) - 0.01168*sin(6*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L);
an.sign                   = (sign(t-an.blackman_delay/an.medium.freq0)-sign(t-an.blackman_delay/an.medium.freq0-an.blackman_L/an.medium.freq0))/2;

an.p_expr                 = 4E4*sin(2*pi*an.medium.freq0*t-6/an.medium.freq0 - an.medium.kappa0 * z)  * an.blwin * an.sign;
an.phi_expr               = int(-an.p_expr/an.medium.rho0*an.medium.freq0,t);

%SIMPLIFIED EXPRESSION VARIABLES
an.phi_sq_expr            = an.phi_expr.^2;
an.lapl_t_expr            = diff(diff(an.phi_sq_expr,t),t)/an.medium.c0^2;
an.nabla_t_sq_expr        = (diff(an.phi_expr,t)/an.medium.c0/an.medium.freq0)^2;
  
%UN-SIMPLIFIED EXPRESSION VARIABLES
an.nabla_x_sq_expr        = (diff(an.phi_expr,z))^2;
an.lapl_x_expr            = diff(diff(an.phi_sq_expr,z),z);

% GET THE VALUES OF EACH VARIABLE FOR EACH TIMEPOINT
an.p                    = double(subs(subs(an.p_expr           , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));
an.phi                  = double(subs(subs(an.phi_expr         , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));
an.phi_sq               = double(subs(subs(an.phi_sq_expr      , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));
an.lapl_t               = double(subs(subs(an.lapl_t_expr      , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));
an.lapl_x               = double(subs(subs(an.lapl_x_expr      , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));
an.nabla_x_sq           = double(subs(subs(an.nabla_x_sq_expr  , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));
an.nabla_t_sq           = double(subs(subs(an.nabla_t_sq_expr  , z, an.domain.zpar(an.depth, sim.title)), t, an.domain.tpar));

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
an.lagrangian       = real(an.nabla_t_sq - an.nabla_x_sq);
an.lagrangian(isnan(an.lagrangian))=0;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.lagrangian, an.lagrangian, an.depth, sim.title);

%% PLOT RANDOM variables
sim.title = ['Lagrangian, X-wave @ 21 deg, Focused'];

an.filename        = [an.file.dirname '/' 'ContrastSrc' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.rnd_var]       = FieldLoad(an.filename);
an.rnd_var          = sim.rnd_var(:,1,1)*0;-(plin.data(:,105,an.depth)).^2/medium.c0^2/medium.rho0/2;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.rnd_var, an.rnd_var, an.depth, sim.title);
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
Y_phi_fourier = abs(fft(sim_var(:,size(sim_var,2)/2,depth),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
f_p = linspace(0,1,Nfft_p/2+1)*F_s;                                                           % [Hz] Frequency sample Number
% Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
S_sim = 20*log10(Y_phi_fourier(Np));

Y_phi_an = abs(fft(an_var,Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
S_an = 20*log10(Y_phi_an(Np));

figure('WindowState','maximized'); hold on; grid on;
title('Frequency spectrum , grad phi')
plot(f_p*1E-6,S_sim-max(S_sim))
plot(f_p*1E-6,S_an-max(S_an))
legend('fourier','an')
xlabel('Frequency [MHz]')
ylabel('Amplitude [dB]')
ylim([-80 0])
title([title_var , ' @  Z = ', num2str(depth), '[mm]'])

set(gca,'FontSize',15)
set(gcf,'Color','white')

OSFactor = 5;
an.domain.tpar_OS = (0:size(an_var,1)*OSFactor-1)*an.domain.dtpar/OSFactor;
figure('WindowState','maximized');hold on; grid on;
title('grad phi computed via different methods')
plot(an.domain.tpar_OS *1E6,                           resampleSINC(sim_var(:,size(sim_var,2)/2,depth)',OSFactor))
plot(an.domain.tpar_OS *1E6, resampleSINC(an_var',OSFactor))
legend('fourier','an')
xlabel('Time [usec]')
ylabel('phi [m^2/sec]')
title([title_var , ' @  Z = ', num2str(depth), '[mm]'])

set(gca,'FontSize',15)
set(gcf,'Color','white')
end

function GENERATE_MAP(an, sim_var, title_var)

%==================================================== UNFOCUSED
figure('WindowState','maximized');
length_z = 1;
imagesc(an.domain.zpar(length_z:an.domain.dims(3) -length_z)*1E3,an.domain.xpar*1E3,squeeze(max(abs(sim_var(:,:,length_z:end-length_z)))));

xlabel('Z [mm]' )
ylabel('X [mm]' )
title(title_var)
colorbar;

file.plot_colour = 'viridis';xAM_Colormaps(file)
set(gca,'FontSize',20)
set(gcf,'Color','white')
%==================================================== FOCUSED
figure('WindowState','maximized')
start_z = 100;
end_z = an.domain.dims(3)-20;
imagesc(an.domain.zpar(start_z:end_z)*1E3,an.domain.xpar*1E3,squeeze(max(abs(sim_var(:,:,start_z:end_z)))));

xlabel('Z [mm]' )
ylabel('X [mm]' )
title([title_var, ', Zoom-In'])
colorbar;
file.plot_colour = 'viridis';xAM_Colormaps(file)

set(gca,'FontSize',20)
set(gcf,'Color','white')
end