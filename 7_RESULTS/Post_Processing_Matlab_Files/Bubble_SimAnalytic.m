% *********************************************************************************
% Simulation of One Scatterer
%
% **********************************************************
% [V_dd_norm] = Single_bubble_sim(p_driv,t_p_driv,c,rho,f)
% **********************************************************
% Computes the second derivative of volume for one scatterer
% by solving the Rayleigh - Plesset Model for microbubbles
%
% Inputs:
% p_driv     real      Array with the values of the driving pressure
% t_p_driv   real      Array with the values of the temporal dimension of p_driv
% c          real      Speed of sound in the medium [m/sec]
% rho        real      Density of medium [kg/m^3]
% f          real      Frequency of the driving pulse [MHz]
%
% Outputs:
% V_dd_norm  real      Second time derivative of volume of the scatterer
%
% ******************
% AM, version 200723
% ******************
%
% *********************************************************************************

function [V_dd_norm] = Single_bubble_sim(p_driv,t_p_driv,c,rho,f)
%% ========================Initialization==============================================================

% -------------Bubble parameters
par.R0 = 2.4e-6 ;                         % [μm],  initial bubble radius R0
par.S_vis = 5.8e-9         ;            % [Pa*sec] , Shell viscosity
par.S_vis = 1.5E-9*exp(8E5*par.R0);

%------------------ Medium parameters (water, Room temperature =20° and 1 atm ambient pressure)
par.P0 = 1.01e+5 ;                             % [Pa], ambient pressure 1 [atm] = 10^5 [Pa]
par.sigma_w  = 0.072 ;                      % [N/m] surface tension
par.sigma_R0 = 0.036;                      % [N/m], surface tension at initial radius R0
par.P_g0 = par.P0+2*par.sigma_R0/par.R0 ;   % [Pa] Initial gas pressure , (Sometimes = P0)
par.gamma = 1.07   ;                           %  (or κ), dimensionless, polytropic exponent of the gas inside the bubble
par.mu = 2e-3 ;                             % [Pa*sec], dynamic viscocity of medium
par.c = c;
par.rho = rho;
par.f = f;
%---------------Marmottant model parameters
par.chi = 0.5   ;                                     % [N/m] shell elasticity
par.R_b = par.R0/sqrt(par.sigma_R0/par.chi+1);     % buckling radius
par.R_r = par.R_b*sqrt(par.sigma_w/par.chi+1);       % upper limit radius after rupture

par.f0 = 1/(2*pi*par.R0*sqrt(par.rho))*sqrt(3*par.gamma*par.P0+ (3*par.gamma-1)*2*par.sigma_R0/par.R0  + 4*par.chi/par.R0);
par.omega = 3*par.P0*par.gamma/par.R0^2/par.rho+4*par.chi/par.R0^3/par.rho;
par.delta = par.R0/par.c + 4*par.mu/(par.rho*par.R0^2*par.omega) + 4*par.S_vis/(par.rho*par.R0^3*par.omega);
par.fresonance = par.f0*sqrt(1-par.delta^2/2)/1E6  %MHz

%% Create p_driv manually
w_driving   = 2*pi*par.f;
gauss_dl    = 12/par.f;
gauss_win   = 6/par.f;
dt = 1/(2*9*f);
t_p_driv = [1:600]*dt;

p_driv   = 1E2*exp( - ((t_p_driv-gauss_dl)/(gauss_win/2)).^2).* sin(w_driving*(t_p_driv-gauss_dl));
%% ========================= Solve ODE with Marmottant ==================================================
T_start = t_p_driv(1) ;
T_end   = t_p_driv(end)   ;

Initial_Cond = [par.R0  ;  0 ]   ;  % initial radius R0 ; initial wall velocity R_dot

% Numerically solve the Rayleigh-Plesset with Marmottant using ode45-solver
% Output: T_m=time, R_m=radius
options  = odeset('RelTol', 1e-12, 'AbsTol', [1e-10 1e-10]', 'Refine', 1);
[T_m, R_m] = ode45 (@(t,R) Bubble_RPAnalytic(t, R, t_p_driv, p_driv, par), [T_start T_end], Initial_Cond, options);

%% ==================== Compare with analytical solution ====================
% Interpolate R_m for even time steps
R = interp1(T_m,R_m,t_p_driv);
% R=R_m;
% p_driv_interp = interp1(t_p_driv,p_driv,T_m);
par.sigma_R = (R(:,1)<par.R_b)*0 + (R(:,1)>par.R_r)*par.sigma_w + (R(:,1)<=par.R_r & R(:,1)>=par.R_b).*par.chi.*(R(:,1).^2./par.R_b.^2-1) ;
P_elas =  2*par.sigma_R./R(:,1) ;
P_vis  =  4*par.S_vis.*R(:,2)./R(:,1).^2 ;
P_gas  = (par.P_g0.*(par.R0./R(:,1)).^(3*par.gamma).*(1-3*par.gamma*R(:,2)/par.c)-par.P0-p_driv'-4*par.mu*R(:,2)./R(:,1) - P_elas-P_vis) ;
R_dd = (P_gas/par.rho-3/2*R(:,2).^2)./R(:,1);
R_d = R(:,2);
R_ = R(:,1);
% R(:,3) = R_dd;

V_dd_norm = (4*pi*R_.*(R_.*R_dd+2.*R_d.^2))';     % Should be divided by the total volume

r = R_;
psc1 = par.rho*R_./r.*((2*R_d.^2 + R_.*R_dd) - R_.^3.*R_d.^2./(2.*r.^3));
psc2 = par.rho*V_dd_norm/4/pi./r';
psc2 = par.rho*V_dd_norm/4/pi./r' - par.rho*R_.d.^2.*(R_/r).^4;
disp(['Pressure @ ' num2str(max(r)) ' [um] is ', num2str(max(psc1*1E-3)) ' [kPa].', num2str(max(psc2*1E-3)), '[kPa].'])
%% ======================= Spectral Response ====================================================
F_s = length(t_p_driv)/(abs(T_end) + abs(T_start));                    % [Hz] , 1/dt
t_interval = [0 T_end]; % time interval to analyse
jj = find(t_p_driv > t_interval(1) & t_p_driv<t_interval(2)); % indices within time interval

% Spectrum of driving pressure
% fft(p_driv(jj))  =        for i = jj-jj(1) ; q(i+1) =  sum(p_driv(jj).*exp(-j*2*pi*(jj-jj(1))*i/length(jj)))  ; end   ;

Nfft_p = 2^nextpow2(2*length(jj));                      % Zero-padding to improve fft performance until length equal to 2N
Y_p = abs(fft(p_driv(jj)./1e6,Nfft_p));                  % Convert the driving pulse to the frequency domain
f_p = F_s/2*linspace(0,1,Nfft_p/2+1);                   % [Hz] Frequency sample Number
% Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 0.5
% Other way to define frequency domain : f_p = Fs*(0:(Nfft_p/2))/Nfft_p
Np = 1:Nfft_p/2+1;                                      % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
P1_p = Y_p(Np);
Sp_p = 20*log10(P1_p);                                  % Convert pressure to amplitude [db]

% Spectrum of radial oscillations
t_data=t_p_driv(jj);
Data=R(jj,1);
Nfft_r = 2^nextpow2(2*length(Data));
Y_r = abs(fft(Data-par.R0,Nfft_r));
f_r = F_s/2*linspace(0,1,Nfft_r/2+1);

Nr = 1:Nfft_r/2+1;
P1_r = Y_r(Nr);
Sp_r = 20*log10(P1_r);
i_ff = find(f_r > f,1)-1;                            % index for fundamental frequency

%% Create plots
fig = figure('WindowState','maximized');
hold on;

% Create plot of pressure driving pulse in respect with time , p(t)
h_f1 = subplot(2,2,1);
h = plot(t_p_driv(jj),p_driv(jj), '-', 'LineWidth', 1.5);
fig_color = get(h,'Color');
xlabel('time [s]'); ylabel('pressure [Pa]'); title('driving pulse'); box on;
hold off;

hold on;
% Create plot of radial oscillation in respect with time , p(t)
h_f2 = subplot(2,2,3);
plot(t_p_driv(jj),R(jj,1), '-', 'LineWidth', 1.5, 'color', fig_color);
xlabel('time [s]'); ylabel('Radius [m]'); title('bubble radius'); box on;
hold off;

% Create plot of frequency spectrum for driving pulse and radial oscillations
h_f3 = subplot(2,2,4);
hold on;
m=plot(f_r,Sp_r-max(Sp_r), '-', 'LineWidth', 1.5, 'color', fig_color);
pr=plot(f_p,Sp_p-max(Sp_p),':','LineWidth', 1, 'color', fig_color);
legend([m, pr], 'Radial', 'Driving pulse' )
ylim([-45 max(Sp_r-max(Sp_r))]);
ylabel('normalized amplitude [dB]'); title('power spectra'); box on;
xlim([0 3*f]);
hold off;

% Create Plot Radius-surface tension curve
h_f4 = subplot(2,2,2);
hold on;

sigma_R = @(R) par.chi*(R.^2/par.R_b^2-1) ;
buck_reg = [min(0.95*par.R_b,min(0.98.*R(:,1))) par.R_b] ;
rupt_reg = [par.R_r max(1.05*par.R_r,max(1.02.*R(:,1)))];
plot(buck_reg, zeros(length(buck_reg),1), ':', 'LineWidth', 1.5, 'color', fig_color); xlabel('radius [m]'); ylabel('surface tension [N/m]')
hold on; plot([par.R_b par.R0 par.R_r], sigma_R([par.R_b par.R0 par.R_r]), ':', 'LineWidth', 1.5, 'color', fig_color);
hold on; plot(rupt_reg, par.sigma_w.*ones(length(rupt_reg),1), ':', 'LineWidth' ,1.5, 'color', fig_color); box on;
xlim([buck_reg(1),rupt_reg(2)])
ylim([0 0.085])

% Label regimes
plot([par.R_b par.R_b], [0 , 0.085],'--', 'color', fig_color);
plot([par.R_r par.R_r], [0 0.085], '--', 'color', fig_color)
text (min(0.95*par.R_b,min(1.01.*R(:,1))), 0.08, 'buckling','FontSize',11 , 'FontWeight', 'bold');
text (max(1.02*par.R_r,max(0.97.*R(:,1))), 0.08, 'rupture','FontSize',11, 'FontWeight', 'bold');
text ((par.R_b+par.R_r)/2, 0.08, 'elastic','FontSize',11, 'FontWeight', 'bold');

% Position of min and max
vert = -0.005 ; hor = - 0.004e-6;
plot(min(R(:,1)),max(sigma_R(min(R(:,1))),0), 'x', 'LineWidth', 1.5, 'MarkerSize', 8, 'color', fig_color);hold on
plot(par.R0,sigma_R(par.R0), 'x', 'LineWidth', 1.5, 'MarkerSize', 8, 'color', fig_color);hold on
plot(max(R(:,1)),min(sigma_R(max(R(:,1))),par.sigma_w), 'x', 'LineWidth', 1.5, 'MarkerSize', 8, 'color', fig_color);hold on
text (min(R(:,1))+hor, max(sigma_R(min(R(:,1))),0)+vert, 'R_{min}' ,'color', fig_color);
text (par.R0+hor, par.sigma_R0+vert, 'R_{0}' , 'color', fig_color);
text (max(R(:,1))+hor, min(sigma_R(max(R(:,1))),par.sigma_w)+vert, 'R_{max}', 'color', fig_color);
title('surface tension');


set(h_f1,'XTickLabel',h_f1.XTick*1e+6) ; h_f1.XLabel.String = 'Time [μsec]';

set(h_f2,'XTickLabel',h_f2.XTick*1e+6) ; h_f2.XLabel.String = 'Time [μsec]';

% set(h_f3,'XTickLabel',h_f3.XTick*1e-6) ; h_f3.XLabel.String = 'Frequency [MHz]';

set(h_f4,'XTickLabel',h_f4.XTick*1e+6) ; h_f4.XLabel.String = 'Radius [μm]';

end

