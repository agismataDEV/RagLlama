
%%
clear an
an.medium.freq0         = 15E6;
an.medium.Fnyq          = domain.Fnyq;
an.medium.rho0          = 1060;
an.medium.c0            = 1480;
an.medium.kappa0        = 2*pi*an.medium.freq0/an.medium.c0;
an.medium.dLambda       = an.medium.c0/an.medium.freq0;

an.domain.start         = domain.start;
an.domain.dtpar         = 1/(2*an.medium.Fnyq*an.medium.freq0);
an.domain.dims          = size(plin.data);
an.domain.tdimpar       = an.domain.dims(1);
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
an.theta       = 10;
%% Symbolics
an.depth                = 258;
an.lat                  = 176;

syms t z x
an.blackman_L             = 6;
an.blackman_delay         = 6;
an.blwin                  = 0.35869 - 0.48829 * cos(2*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) + 0.14128*cos(4*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) - 0.01168*sin(6*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L);
an.sign                   = (sign(t-an.blackman_delay/an.medium.freq0)-sign(t-an.blackman_delay/an.medium.freq0-an.blackman_L/an.medium.freq0))/2;

% an.p_expr                 = 4E4*sin(2*pi*an.medium.freq0*t-6/an.medium.freq0 - an.medium.kappa0 * z)  * an.blwin * an.sign;

an.gauss_dl                  = 6.0/an.medium.freq0;
an.gauss_win                 = 3.0/an.medium.freq0;
w                            = 2 * pi * an.medium.freq0;

A                         = ((t-an.gauss_dl)/(an.gauss_win/2));
% an.p_expr                 = 4E5*exp( - A.^2).* sin(w*(t-an.gauss_dl) - an.medium.kappa0*sind(an.theta)*x - an.medium.kappa0*cosd(an.theta)*z);
% an.p_expr                 = 1E5*exp( - A.^2).*(sin(w*(t-gauss_dl) - an.medium.kappa0*sind(theta)*x + an.medium.kappa0*cosd(theta)*z)*-A/(gauss_win/2*an.medium.freq0) + cos(w*(t-gauss_dl)- an.medium.kappa0*sind(theta)*x + an.medium.kappa0*cosd(theta)*z));
an.p_expr                   = 4E5*(sin(w*(t-an.gauss_dl) - an.medium.kappa0*sind(an.theta)*x - an.medium.kappa0*cosd(an.theta)*z) + sin(w*(t-an.gauss_dl) + an.medium.kappa0*sind(an.theta)*x- an.medium.kappa0*cosd(an.theta)*z));
% an.p_expr                   = 3.3E5*(sin(w*(t-gauss_dl) - an.medium.kappa0*sind(theta)*x - an.medium.kappa0*cosd(theta)*z));
an.phi_expr                 = int(-an.p_expr,t)/an.medium.rho0*an.medium.freq0;

%UN-SIMPLIFIED EXPRESSION VARIABLES
an.nabla_z_sq_expr        = (diff(an.phi_expr,z))^2/an.medium.freq0^2;
an.nabla_x_sq_expr        = (diff(an.phi_expr,x))^2/an.medium.freq0^2;
an.nabla_t_sq_expr        = (an.p_expr/an.medium.c0/an.medium.rho0)^2;%(diff(an.phi_expr,t)/an.medium.c0/an.medium.freq0)^2;
an.lagr_expr              =  an.medium.rho0/2 * ( an.nabla_z_sq_expr + an.nabla_x_sq_expr - an.nabla_t_sq_expr);

%UN-SIMPLIFIED EXPRESSION SOURCE TERM
an.lagr_first_xyz_deriv   = diff(an.lagr_expr,z) + diff(an.lagr_expr,x);
an.lagr_z_expr            = diff(an.lagr_first_xyz_deriv,z);
an.lagr_x_expr            = diff(an.lagr_first_xyz_deriv,x);
an.lagr_t_expr            = diff(diff(an.lagr_expr,t),t)/an.medium.c0^2;%(diff(an.phi_expr,t)/an.medium.c0/an.medium.freq0)^2;
an.source_term_expr       = an.lagr_z_expr + an.lagr_x_expr + an.lagr_t_expr;

an.lagr_pres_x_expr       = diff(diff(an.phi_expr^2,x),x);
an.lagr_pres_z_expr       = diff(diff(an.phi_expr^2,z),z);
an.lagr_pres_t_expr       = diff(diff(an.phi_expr^2,t),t)/an.medium.c0^2;
an.lagr_pres_expr         = an.medium.rho0/4*(an.lagr_pres_z_expr + an.lagr_pres_x_expr + an.lagr_pres_t_expr)/an.medium.freq0^2;

%SIMPLIFIED EXPRESSION VARIABLES
an.phi_sq_expr_simpl                = an.phi_expr.^2;
an.nabla_t_sq_expr_simpl            = diff(diff(an.phi_sq_expr_simpl,t),t)/an.medium.c0^2;
an.nabla_x_sq_expr_simpl            = diff(diff(an.phi_sq_expr_simpl,z),z);

an.plin = plin;
an.pnl  = pnl;
an.dslice = dslice;
an.file   = file;

%%
% GET THE VALUES OF EACH VARIABLE FOR EACH TIMEPOINT
an.p                    = double(subs(subs(an.p_expr           , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.phi                  = double(subs(subs(an.phi_expr         , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.phi_sq               = double(subs(subs(an.phi_sq_expr_simpl, [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.nabla_x_sq           = double(subs(subs(an.nabla_x_sq_expr  , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.nabla_z_sq           = double(subs(subs(an.nabla_z_sq_expr  , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.nabla_t_sq           = double(subs(subs(an.nabla_t_sq_expr  , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.lagrangian           = double(subs(subs(an.lagr_expr        , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.source_term          = double(subs(subs(an.source_term_expr , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
an.lagr_pres            = double(subs(subs(an.lagr_pres_expr   , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
%% Implement the window and the phase shift based on the element of intererst

%% Lagrangian
an.lagrangian           = double(subs(subs(an.lagr_expr        , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
[var,an] = Compare_Signals(an, 13, an.lagrangian, 0); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
an.lagrangian = var;

%% Source Term
an.source_term          = double(subs(subs(an.source_term_expr , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
[var,an] = Compare_Signals(an, 13, an.source_term, 0); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
an.source_term = var;

%% Pressure Term
an.lagr_pres            = double(subs(subs(an.lagr_pres_expr   , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
[var,an] = Compare_Signals(an, 13, an.lagr_pres, 0); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
an.lagr_pres = var;
%% Pressure
an.p                    = double(subs(subs(an.p_expr           , [z x], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), t, an.domain.tpar));
[var,an] = Compare_Signals(an, 13, an.p, 0); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
an.p = var;

%%  ====================================================== P ======================================================
sim.title       = ['Incident Pressure Field @ 21$^{\circ}$, Focused'];
sim.ylabel_var  = ['Pressure [kPa]'];
an.filename     = [an.file.dirname '/' 'TESTNeumann' int2string_ICS(0) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
% [sim.p]         = FieldLoad(an.filename);
an.OSFactor     = 4;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.p*1E-3, an.p'*1E-3, an.depth, sim.title, sim.ylabel_var);
GENERATE_MAP(an, sim.p*1E-3, sim.title, 'Pressure [kPa]');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====================================================== Lagrangian ======================================================
close all;
sim.title           = ['Lagrangian Density @ ', num2str(an.theta),'$^{\circ}$, Focused'];
sim.ylabel_var      = ['Lagrangian Density [kg/m$\cdot$s$^2$]'];
an.filename         = [an.file.dirname '/' 'Lagrangian' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
% [sim.lagrangian]    = FieldLoad(an.filename);
an.OSFactor         = 1;

[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.lagrangian, an.lagrangian', an.depth, sim.title,sim.ylabel_var);
GENERATE_MAP(an, sim.lagrangian, sim.title,sim.ylabel_var);

%% ====================================================== Source Term due to Lagrangian ======================================================
% close all;
sim.title                                                       = ['Source Term due to Lagrangian Density @ ', num2str(an.theta),'$^{\circ}$'];
sim.ylabel_var                                                  = ['Amplitude normalized with $\lambda^2$ [Pa]'];
an.filename                                                     = [an.file.dirname '/' 'LagrangianPressure' int2string_ICS(2) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.lagrangian_source_term]                                    = FieldLoad(an.filename);
an.lagrangian_source_term                                       = real(double(an.source_term));
sim.Lagr_ST                                                     = sim.lagrangian_source_term/an.medium.dLambda ^2; % Normalization Factor
an.OSFactor                                                     = 4;
[an.FreqRange,sim.FreqSpect,an.FreqSpect]                       = Frequency(an, sim.Lagr_ST, an.lagrangian_source_term', an.depth, sim.title, sim.ylabel_var);
                                                                  GENERATE_MAP(an, sim.Lagr_ST, sim.title, sim.ylabel_var);

%% ====================================================== Source Term due to MediumNL ======================================================
% close all;
sim.title                                                       = ['Source Term due to medium nonlinearities @ ', num2str(an.theta),'$^{\circ}$'];
sim.ylabel_var                                                  = ['Amplitude normalized with $\lambda^2$ [Pa]'];
an.filename                                                     = [an.file.dirname '/' 'MediumNLCST' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.mediumnl_source_term]                                      = FieldLoad(an.filename);
an.mediumnl_source_term                                         = real(double(an.source_term));
sim.MNL_ST                                                      = sim.mediumnl_source_term/an.medium.dLambda ^2; % Normalization Factor
an.OSFactor                                                     = 4;
[an.FreqRange,sim.FreqSpect,an.FreqSpect]                       = Frequency(an, sim.MNL_ST, an.mediumnl_source_term', an.depth, sim.title, sim.ylabel_var);
                                                                  GENERATE_MAP(an, sim.MNL_ST, sim.title, sim.ylabel_var);

%% PLOT RANDOM variables of interest
sim.title = ['Pressure due to Lagrangian @ 21$^{\circ}$'];
sim.ylabel_var                                                  = ['Pressure due to Lagrangian [Pa]'];

an.filename        = [an.file.dirname '/' 'ContrastSrc' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
% [sim.rnd_var]       = FieldLoad(an.filename);
an.rnd_var          = real(double(an.lagr_pres));
an.OSFactor                                                     = 4;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.rnd_var, an.rnd_var', an.depth, sim.title, sim.ylabel_var);
GENERATE_MAP(an, sim.rnd_var, sim.title, sim.ylabel_var);

%%===================================================================================================================================================
%% ------------------------------------------------------ GENERAL FUNCTIONS -------------------------------------------------------------------------
%%===================================================================================================================================================
%%
function [fourier] = FieldLoad(filename)

output_nonl         = load_ICS_slice(filename);
fourier             = squeeze(output_nonl.data);

end
%%
function  [f_p,S_sim,S_an] = Frequency(an, sim_var, an_var ,depth, title_var, ylabel_var)

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

an.domain.tpar_OS = (0:size(an_var,1)*an.OSFactor-1)*an.domain.dtpar/an.OSFactor;

figure('WindowState','maximized');hold on; grid on;
plot(an.domain.tpar_OS *1E6, resampleSINC(sim_var(:,an.lat,depth)',an.OSFactor),'--','Color','blue','LineWidth',2)
plot(an.domain.tpar_OS *1E6, resampleSINC(an_var',an.OSFactor),'k-','LineWidth',2)

legend('INCS','Analytical','Interpreter','latex')
xlabel('Time [$\mu$s]','Interpreter','latex')
ylabel(ylabel_var,'Interpreter','latex')
title([title_var , ' @ (X,Y,Z) =(', num2str(an.domain.xpar(an.lat)*1E3),', ', num2str(0),', ', num2str(an.domain.zpar(an.depth)*1E3),') [mm]'],'Interpreter','latex')

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',15)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)
end
%%
function GENERATE_MAP(an, sim_var, title_var, ylabel_val)

%==================================================== UNFOCUSED
figure('WindowState','maximized');
length_z = 1;
an.domain.zpar_OS = linspace(an.domain.zpar(1),an.domain.zpar(end),an.OSFactor*(length(an.domain.zpar))+1);
an.domain.xpar_OS = linspace(an.domain.xpar(1),an.domain.xpar(end),an.OSFactor*(length(an.domain.xpar))+1);

[Z,X] = meshgrid(an.domain.zpar,an.domain.xpar);
plot_var = squeeze(max(abs(sim_var)));
[Zq,Xq] = meshgrid(an.domain.zpar_OS ,an.domain.xpar_OS);

plot_varINTERP = interp2(Z,X,plot_var,Zq,Xq,'spline');

imagesc(an.domain.zpar_OS*1E3,an.domain.xpar_OS*1E3,plot_varINTERP);

xlabel('Z [mm]','Interpreter','latex')
ylabel('X [mm]','Interpreter','latex' )
title(title_var,'Interpreter','latex')
c=colorbar;
c.Label.String = ylabel_val;
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
start_z = 300;
end_z = length(an.domain.zpar_OS)-300;
imagesc(an.domain.zpar_OS(start_z:end_z)*1E3,an.domain.xpar_OS*1E3,plot_varINTERP(:,start_z:end_z));

xlabel('Z [mm]','Interpreter','latex' )
ylabel('X [mm]','Interpreter','latex' )
title([title_var, ', Zoom-In'],'Interpreter','latex')

c=colorbar;
c.Label.String = ylabel_val;
c.Label.Interpreter='latex';
c.TickLabelInterpreter='latex';
file.plot_colour = 'viridis';xAM_Colormaps(file)

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',20)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)
end
%%
function plot_time_frame(an,frame)
%%
i_end = an.domain.tdimpar;
plot_value = (20*log10(abs(an.pnl.xwave_data)));
dBVALUES = 0.3*max(max(max(plot_value)));

%% ========================= Initialize first frame=========================
% This is done to speedup the process of making images , the most time is lost in making the plots
num_frames = i_end;
F(num_frames) = struct('cdata',[],'colormap',[]);

f2 = figure('WindowState','maximized');
ax = gca;
p=imagesc(ax,an.domain.zpar,an.domain.xpar,squeeze(plot_value(1,:,:)) );
hold on;
caxis(ax,[max(max(max(plot_value)))-dBVALUES max(max(max(plot_value)))])
title_txt = 'Left-Aperture Pressure Field, $t_i$ = ';
title(ax,[title_txt,num2str(1)],'Interpreter','latex')
xlabel(ax,an.dslice.xlabel,'Interpreter','latex')
ylabel(ax,an.dslice.ylabel,'Interpreter','latex')

xlim(ax,[min(an.domain.zpar) max(an.domain.zpar)])
ylim(ax,[min(an.domain.xpar) max(an.domain.xpar)])

xAM_Colormaps(an.file)
c = colorbar;
c.Label.String = 'Pressure [dB]';
c.Label.Interpreter ='latex';
c.TickLabelInterpreter = 'latex';

set(ax,'YDir','normal')
if an.dslice.savedim(an.dslice.num) =='z'; set(ax,'XDir','reverse');set(ax,'YDir','reverse') ;  camorbit(90,180);    camroll(180) ; end
axis image;
set(gcf, 'Color', 'w');
set(gca,'FontSize',18)
ax.TickLabelInterpreter = 'latex';
set(findall(gca,'-property','FontSize'),'FontSize',30)

for i=frame:frame % Campure a specific time frame and plot it in a figure
    set(ax.Title,'String',[title_txt,num2str(i)],'Interpreter','latex')
    set(p,'XData',an.domain.zpar,'YData',an.domain.xpar,'CData',squeeze(plot_value(i,:,:)))
    drawnow
end
plot([an.domain.zpar(1) an.domain.zpar(end)],[0 0],'--','Color','white','LineWidth',2)
%Capture the indexes of the max values in this plot
[max_val,max_idx] = max(squeeze(plot_value(i,:,:))');

% Plot the values in a figure
f3 = figure('WindowState','maximized');
plot(an.domain.zpar(max_idx),an.domain.xpar,'LineWidth',2)

xlabel(an.dslice.xlabel,'Interpreter','latex')
ylabel(an.dslice.ylabel,'Interpreter','latex')
axis image;
xlim([an.domain.zpar(1) an.domain.zpar(end)])

% Add 2 points in the same figure
hold on;
P_ZX = [an.domain.zpar(max_idx) ; an.domain.xpar]';
P1 = P_ZX(160,:);
P2 = P_ZX(191,:);
plot(P1(:,1),P1(:,2),'ro','LineWidth',2','MarkerSize',5)
plot(P2(:,1),P2(:,2),'ro','LineWidth',2','MarkerSize',5)

%Find many differences
P1_many = P_ZX(160,:);
P2_many = P_ZX(191,:);
Angle_many=90-atand((P2_many(:,2)-P1_many(:,2))./(P2_many(:,1)-P1_many(:,1)));
title(['Max pressure location in X-Z plane, ' num2str(an.theta) ' deg, @ Y = 0 mm '],'Interpreter','latex')
%     title(['Average angle between the 2 points, ',num2str(mean(Angle_many)) ,' deg'])
legend(['Angle between 2 points: ',num2str(mean(Angle_many)) ,' deg'],'Interpreter','latex')
plot([an.domain.zpar(1) an.domain.zpar(end)],[0 0],'--','Color','black','LineWidth',2)

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(findall(gca,'-property','FontSize'),'FontSize',40)

set(gcf, 'Color', 'w');
set(gca,'FontSize',18)
end
%%
function [an_press,an] = Compare_Signals(an, element_n, var, plot_val)
%%
[x_el,arr,td_x,apod_x,apod_right,apod_left] = xAM_Generate_TD_APOD(an.theta);
dist_cross_x_axis = an.domain.zpar(an.depth);
if (an.theta>0)
    element_x = abs(x_el(element_n))*1E-3; % m
    dist_cross_x_axis = element_x*tand(90-an.theta); % m
    [val,an.depth]= min(abs(an.domain.zpar-dist_cross_x_axis));
else
    td_x = zeros(32,1);element_x =0;
end
travel_dist_xwave = dist_cross_x_axis*cosd(an.theta) + element_x*sind(an.theta) + an.domain.start(3)*an.domain.dxpar; % m
an_td = travel_dist_xwave/an.medium.c0 + td_x(32 - element_n +1)*1E-6; % There is some small difference because the dist that cross x is not always a gridpoint in the computational domain
% integer_td = round(an_td/an.domain.dtpar);
% an_td = integer_td*an.domain.dtpar;
frame = round((an_td + an.gauss_dl)/an.domain.dtpar );
plot_time_frame(an,frame)

window = an.gauss_dl * an.medium.freq0 * 2 * an.medium.Fnyq + an.depth/an.medium.c0;
len_window  = (an.gauss_win + an.gauss_dl)  * an.medium.freq0 * 2 * an.medium.Fnyq;

an.rect_win  = [zeros(round(an_td/an.domain.dtpar),1); rectwin(len_window); zeros(round(an.domain.tdimpar-an_td/an.domain.dtpar - len_window),1)];
an_press  = real(ifft(fft(var).*exp(-1i*an.omega*an_td))).* an.rect_win';

if (plot_val)
    figure('WindowState','maximized');
    hold on;plot(an.domain.xpar*1E3,max(abs(an.plin.data(:,:,an.depth)))*1E-3,'LineWidth',2)

    ax = gca;
    xlabel('X [mm]','Interpreter','latex')
    ylabel('Pressure [kPa]','Interpreter','latex')
    title(['Lateral profile @ Z =', num2str(an.domain.zpar(an.depth)*1E3),' [mm]'],'Interpreter','latex')

    ax.TickLabelInterpreter = 'latex';
    set(findall(gca,'-property','FontSize'),'FontSize',40)
    set(gcf,'Color','white')
else
    close all
end
end