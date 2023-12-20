function xAM_Compare_Analytical()
%% Set initial parameters for comparison
clear an
an.plin = plin;
an.pnl  = pnl;
an.dslice = dslice;
an.file   = file;

an.medium.freq0         = medium.freq0;
an.medium.Fnyq          = domain.Fnyq;
an.medium.Fs            = domain.Fnyq*medium.freq0*2;
an.medium.rho0          = 1060;
an.medium.c0            = 1480;
w                       = 2 * pi * an.medium.freq0;
an.medium.kappa0        = w/an.medium.c0;
an.medium.dLambda       = an.medium.c0/an.medium.freq0;

an.OSFactor             = 1;
an.domain.start         = domain.start;
an.domain.dtpar         = 1/(2*an.medium.Fnyq*an.medium.freq0*an.OSFactor);
an.domain.dims          = size(plin.data);
an.domain.tdimpar       = an.domain.dims(1)*an.OSFactor;
an.domain.tpar          = [0:an.domain.tdimpar-1]*an.domain.dtpar ;
an.domain.tdimpar       = length(an.domain.tpar);

an.domain.dxpar         = an.domain.dtpar*an.medium.c0;
an.domain.zpar          = domain.par{3}*1E-3;
an.domain.xstart        =  -an.domain.dims(2)/2;
an.domain.xpar          = (an.domain.xstart+(0:an.domain.dims(2)-1))*an.domain.dxpar;%an.domain.xpar = domain.par{1};
an.domain.ypar          = domain.par{2}*1e-3;

an.domain.beam          = 0;
an.file.focalplanename  = [ 'y' int2string_ICS(2) ];

Fs          = 1/an.domain.dtpar;  dF = Fs/an.domain.tdimpar;
if mod(an.domain.tdimpar,2) == 0 % Nt is even
    f = -Fs/2 : dF : Fs/2-dF;
else            % Nt is odd
    f = [sort(-1*(dF:dF:Fs/2)) (0:dF:Fs/2)];
end
an.omega       = 2*pi*f;
an.omega       = [an.omega(an.omega>=0) an.omega(an.omega<0)];
an.theta       = 20;

an.depth                = 276;
an.lat                  = round(size(plin.data,2)/2+1);
%% Derive the analytical expressions in symbolic format
% Everything is normalized with center frequency
syms t tnorm z znorm x xnorm
an.blackman_L             = 12;
an.blackman_delay         = 6;
an.blwin                  = 0.35869 - 0.48829 * cos(2*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) + 0.14128*cos(4*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L) - 0.01168*sin(6*pi*(an.medium.freq0*t-an.blackman_delay)/an.blackman_L);
an.sign                   = (sign(t-an.blackman_delay/an.medium.freq0)-sign(t-an.blackman_delay/an.medium.freq0-an.blackman_L/an.medium.freq0))/2;

% an.p_expr              = 4E4*sin(2*pi*an.medium.freq0*t-6/an.medium.freq0 - an.medium.kappa0 * z)  * an.blwin * an.sign;

an.gauss_dl                  = 10/an.medium.freq0 * an.OSFactor;
an.gauss_win                 = 11/an.medium.freq0* an.OSFactor;

t                         = tnorm;
x                         = xnorm;
z                         = znorm;

% This holds only in the middle where the 2
% planewaves arrive at the same time.
an.power =   2 ;
A                         = ((t-an.gauss_dl)/(an.gauss_win/2));
an.psin                   = sin(w*(t-an.gauss_dl) - an.medium.kappa0*sind(an.theta)*x - an.medium.kappa0*cosd(an.theta)*z)+ sin(w*(t-an.gauss_dl) + an.medium.kappa0*sind(an.theta)*x- an.medium.kappa0*cosd(an.theta)*z);
an.p_expr                   = 0.85E5*an.psin;
an.phi_expr               = int(-an.p_expr,t)/an.medium.rho0;

%UN-SIMPLIFIED EXPRESSION VARIABLES
an.nabla_z_sq_expr        = (diff(an.phi_expr,znorm)).^2;
an.nabla_x_sq_expr        = (diff(an.phi_expr,xnorm)).^2;
an.nabla_t_sq_expr        = (diff(an.phi_expr,tnorm)/an.medium.c0).^2;
% an.nabla_t_sq_expr         = (an.p_expr/an.medium.c0/an.medium.rho0).^2;
an.lagr_expr              =  an.medium.rho0/2 * ( an.nabla_z_sq_expr + an.nabla_x_sq_expr -  an.nabla_t_sq_expr);

%SIMPLIFIED EXPRESSION VARIABLES
% an.phi_sq_expr_simpl                = an.phi_expr.^2;
% an.nabla_x_sq_expr_simpl            = diff(diff(an.phi_sq_expr_simpl,xnorm),xnorm)/an.medium.freq0^2;
% an.nabla_z_sq_expr_simpl            = diff(diff(an.phi_sq_expr_simpl,znorm),znorm)/an.medium.freq0^2;
% an.nabla_t_sq_expr_simpl            = diff(diff(an.phi_sq_expr_simpl,tnorm),tnorm)/an.medium.c0^2/an.medium.freq0^2;
% an.lagr_expr                        = an.medium.rho0/4 * ( an.nabla_z_sq_expr_simpl + an.nabla_x_sq_expr_simpl - an.nabla_t_sq_expr_simpl);

%SOURCE TERM EXPRESSION
an.lagr_z_expr            = diff(diff(an.lagr_expr,znorm),znorm);
an.lagr_x_expr            = diff(diff(an.lagr_expr,xnorm),xnorm);
an.lagr_t_expr            = diff(diff(an.lagr_expr,tnorm),tnorm)/an.medium.c0^2;      %(diff(an.phi_expr,t)/an.medium.c0/an.medium.freq0)^2;
an.source_term_expr       = (an.lagr_z_expr + an.lagr_x_expr + an.lagr_t_expr);

% PRESSURE ESTIMATE BY HAMILTON
an.lagr_pres_x_expr       = diff(diff(an.phi_expr^2,xnorm),xnorm)/an.medium.freq0^2;
an.lagr_pres_z_expr       = diff(diff(an.phi_expr^2,znorm),znorm)/an.medium.freq0^2;
an.lagr_pres_t_expr       = diff(diff(an.phi_expr^2,tnorm),tnorm)/an.medium.c0^2/an.medium.freq0^2;
an.lagr_pres_expr         = an.medium.rho0/4*( an.lagr_pres_t_expr + an.lagr_pres_x_expr + an.lagr_pres_z_expr);
%
% % WINDOW FOR SIN FUNCTION
an.len_before = floor(an.gauss_dl * an.medium.freq0 * 2 * an.medium.Fnyq);
an.len_window =  floor(an.gauss_win * an.medium.freq0 * 2 * an.medium.Fnyq );
an.len_after  = length(an.domain.tpar)  - an.len_window - an.len_before;
an.rect_win  = [zeros(an.len_before,1); rectwin(an.len_window); zeros(an.len_after,1)];

%% Source Term for lagrangian density
an.OSFactor = 1;
an.td = 0.18E-8;
an.source_term          = double(subs(subs(an.source_term_expr , [znorm xnorm], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), tnorm, an.domain.tpar));
an.source_term          = an.source_term.* an.rect_win';
[var,an] = Compare_Signals(an, an.source_term, 0); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
an.source_term = var;

%% Pressure Term due to local nonlinearities
an.OSFactor = 1;
an.td = 0.18E-8;
an.lagr_pres            = double(subs(subs(an.lagr_pres_expr   , [znorm xnorm], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), tnorm, an.domain.tpar));
an.lagr_pres            = an.lagr_pres.* an.rect_win';
[var,an] = Compare_Signals(an, an.lagr_pres, 0); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
% sim.Lagr_ST=pnl.xwave_data-plin.data;
% [an.FreqRange,sim.FreqSpect,an.FreqSpect]                       = Frequency(an, sim.Lagr_ST, an.lagr_pres', an.depth, '', '');
an.lagr_pres = var;
%% Incident Pressure Field
an.OSFactor = 1;
an.td = 0.18E-8;
an.p                    = double(subs(subs(an.p_expr           , [znorm xnorm], [an.domain.zpar(an.depth), an.domain.xpar(an.lat)]), tnorm, an.domain.tpar));
an.p                    = an.p.* an.rect_win';
[var,an] = Compare_Signals(an, an.p, 1); % Include all variables from an structure, element number, variable of interert, plot (1)/no-plot(0)
an.p = var;

%%  ====================================================== P ======================================================
sim.title       = ['Incident Pressure Field, ',num2str(an.theta),'$^{\circ}$'];
sim.ylabel_var  = ['Normalized Pressure [dB]'];
an.filename     = [an.file.dirname '/' 'TESTNeumann' int2string_ICS(0) '_y_002' int2string_ICS(an.domain.beam)];
% [sim.p]         = FieldLoad(an.filename);
sim.p(sim.p==0)                                     = 1E-30;
an.OSFactor     = 10;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.p, an.p', an.depth, sim.title, sim.ylabel_var);
sim_var = squeeze((20*log10(max(abs(hilbert(sim.p))))));
GENERATE_MAP(an, sim_var , sim.title,sim.ylabel_var);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ====================================================== Lagrangian ======================================================
% close all;
sim.title           = ['Lagrangian Density @ ', num2str(an.theta),'$^{\circ}$, Focused'];
sim.ylabel_var      = ['Lagrangian Density [kg/m$\cdot$s$^2$]'];
an.filename         = [an.file.dirname '/' 'Lagrangian' int2string_ICS(1) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
[sim.lagrangian]    = FieldLoad(an.filename);
sim.lagrangian = sim.Lagr_ST;
sim.lagrangian(sim.lagrangian==0)                                     = 1E-30;
an.OSFactor         = 4;

[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an, sim.lagrangian, an.lagrangian', an.depth, sim.title,sim.ylabel_var);
GENERATE_MAP(an, 20*log10(abs(sim.lagrangian)), sim.title,sim.ylabel_var);

%% ====================================================== Source Term due to Lagrangian ======================================================
% close all;
sim.title                                                       = ['Contrast Source Term, ', num2str(an.theta),'$^{\circ}$'];
sim.ylabel_var                                                  = ['Amplitude [dB/mm$^2$]'];
an.filename                                                     = [an.file.dirname '/' 'LagrangianPressure' int2string_ICS(4) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
% [sim.lagrangian_source_term]                                    = FieldLoad(an.filename);
A=((an.domain.tpar-5.2/an.medium.freq0)/(9.8/an.medium.freq0/2));
q = exp(-A.^12);
an.lagrangian_source_term                                       = real(double(an.source_term)).*circshift(q,360);
sim.Lagr_ST                                                     = sim.lagrangian_source_term(1:end,:,:)/(an.medium.dLambda)^2*1E-11 ; % Normalization Factor
sim.Lagr_ST(sim.Lagr_ST==0)                                     = 1E-30;
an.OSFactor                                                     = 10;
[an.FreqRange,sim.FreqSpect,an.FreqSpect]                       = Frequency(an, sim.Lagr_ST, an.lagrangian_source_term', an.depth, sim.title, sim.ylabel_var,[1.5 2.5]*an.medium.freq0);
%%
sim_var  = 20*log10(maxhilbert_vec(sim.Lagr_ST,an.domain.dtpar,1.5*an.medium.freq0,2.5*an.medium.freq0,4));
%%
sim_var                                                         = squeeze(20*log10(max(abs(hilbert(sim.Lagr_ST)))));
GENERATE_MAP(an, sim_var, sim.title, sim.ylabel_var);

%% ====================================================== Source Term due to MediumNL ======================================================
% close all;
sim.title                                                       = ['Source Term due to medium nonlinearities @ ', num2str(an.theta),'$^{\circ}$'];
sim.ylabel_var                                                  = ['Amplitude [dB/mm$^2$]'];
an.filename                                                     = [an.file.dirname '/' 'MediumNLCST' int2string_ICS(4) '_' an.file.focalplanename int2string_ICS(an.domain.beam)];
% [sim.mediumnl_source_term]                                      = FieldLoad(an.filename);
an.mediumnl_source_term                                         = real(double(an.source_term));
sim.MNL_ST                                                      = sim.mediumnl_source_term/(an.medium.dLambda)^2*1E-11; % Normalization Factor
sim.MNL_ST(sim.MNL_ST==0)                                        = 1E-30;
an.OSFactor                                                     = 10;
[an.FreqRange,sim.FreqSpect,an.FreqSpect]                       = Frequency(an, sim.MNL_ST, an.mediumnl_source_term', an.depth, sim.title, sim.ylabel_var,[1.5 2.5]*an.medium.freq0);
%%
sim_var                                                         = squeeze(max(20*log10(abs(hilbert(sim.MNL_ST)))));
GENERATE_MAP(an, sim_var, sim.title, sim.ylabel_var);

%% PLOT RANDOM variables of interest
sim.title = ['Pressure Field $\vert$ Local Nonlinearities, ',num2str(an.theta),'$^{\circ}$'];
sim.ylabel_var                                                  = ['Normalized Pressure [dB]'];

an.filename        = [an.file.dirname '/' 'ContrastSrc' int2string_ICS(4) '_y_002' int2string_ICS(an.domain.beam)];
% [sim.rnd_var]       = FieldLoad(an.filename);
[sim.rnd_var_db]       = sim.rnd_var;
sim.rnd_var_db(sim.rnd_var_db==0)                                     = 1E-30;
an.rnd_var          = real(double(an.lagr_pres))*an.medium.freq0^2;
an.OSFactor         = 10;
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an,sim.rnd_var_db , -an.rnd_var', an.depth, sim.title, sim.ylabel_var,[1.5*an.medium.freq0 2.5*an.medium.freq0]);
%%
sim_var  = 20*log10(maxhilbert_vec(sim.rnd_var_db,an.domain.dtpar,0.5*an.medium.freq0,3.1*an.medium.freq0,4));
%%
sim_var  = squeeze(20*log10(max(abs(hilbert(sim.rnd_var_db)))));
%%
GENERATE_MAP(an, sim_var, sim.title, sim.ylabel_var);

%%
[an.FreqRange,sim.FreqSpect,an.FreqSpect] = Frequency(an,sim.rnd_var_db , -an.rnd_var', an.depth, sim.title, 'test',[medium.freq0 3*medium.freq0]);
end

%% ------------------------------------------------------ GENERAL FUNCTIONS -------------------------------------------------------------------------
%%===================================================================================================================================================
%%
function [fourier] = FieldLoad(filename)

output_nonl         = load_ICS_slice(filename);
fourier             = squeeze(output_nonl.data);

end
%%
function  [f_p,S_sim,S_an] = Frequency(an, sim_var, an_var ,depth, title_var, ylabel_var, BP_lim)
%%
t_start = 1;
t_end = length(an.domain.tpar);
F_s = 1/(2*an.domain.dtpar);
% Calculate the frequency spectrum based on a sampling frequency and the pressure values
Nfft_p = 2^nextpow2(2*length(an.domain.tpar(t_start:t_end)));                      % Zero-padding to improve fft performance until length equal to 2N
Y_phi_fourier = abs(fft(sim_var(t_start:t_end,an.lat,depth),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
f_p = linspace(0,1,Nfft_p/2+1)*F_s;                                                           % [Hz] Frequency sample Number
% Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
S_sim = 20*log10(Y_phi_fourier(Np));

Y_phi_an = abs(fft(an_var(t_start:t_end),Nfft_p))/Nfft_p;                              % Convert the driving pulse to the frequency domain
S_an = 20*log10(Y_phi_an(Np));

figure('WindowState','maximized'); hold on; grid on;
plot(f_p*1E-6,S_an,'k-','LineWidth',2)
plot(f_p*1E-6,S_sim,'--','Color','blue','LineWidth',2)
legend('Analytical','INCS','Interpreter','latex')
xlabel('Frequency [MHz]','Interpreter','latex')
ylabel('Amplitude [dB]','Interpreter','latex')
ylim([-40 40])
title([title_var , ' @ (X,Y,Z) =(', num2str(an.domain.xpar(an.lat)*1E3),', ', num2str(0),', ', num2str(an.domain.zpar(an.depth)*1E3),') [mm]'],'Interpreter','latex')
%%

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(gca,'FontSize',15)
set(gcf,'Color','white')
set(findall(gca,'-property','FontSize'),'FontSize',30)

an.domain.tpar_OS = (0:size(an_var,1)*an.OSFactor-1)*an.domain.dtpar/an.OSFactor;
sim.domain.tpar_OS = (0:size(sim_var,1)*an.OSFactor-1)*an.domain.dtpar/an.OSFactor;

figure('WindowState','maximized');hold on; grid on;
if (nargin<7)
    plot(an.domain.tpar_OS *1E6, resampleSINC((an_var'),an.OSFactor),'k-','LineWidth',2)
    plot(sim.domain.tpar_OS *1E6, resampleSINC((sim_var(:,an.lat,depth)'),an.OSFactor),'--','Color','blue','LineWidth',2)
else
    plot(an.domain.tpar_OS *1E6, interpft(bandpass(an_var',BP_lim,an.medium.Fs),an.OSFactor*length(an.domain.tpar)),'k-','LineWidth',2)
    plot(sim.domain.tpar_OS *1E6, interpft(bandpass(sim_var(:,an.lat,depth)',BP_lim,an.medium.Fs),an.OSFactor*length(an.domain.tpar)),'--','Color','blue','LineWidth',2)
end

legend('Analytical','INCS','Interpreter','latex')
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

%% ==================================================== UNFOCUSED
figure('WindowState','maximized');
an.domain.zpar_OS = linspace(an.domain.zpar(1),an.domain.zpar(end),an.OSFactor*(length(an.domain.zpar))+1);
an.domain.xpar_OS = linspace(an.domain.xpar(1),an.domain.xpar(end),an.OSFactor*(length(an.domain.xpar))+1);

[Z,X] = meshgrid(an.domain.zpar,an.domain.xpar);
plot_var =sim_var;
[Zq,Xq] = meshgrid(an.domain.zpar_OS ,an.domain.xpar_OS);

plot_varINTERP = interp2(Z,X,plot_var,Zq,Xq,'spline');

imagesc(an.domain.xpar_OS*1E3,an.domain.zpar_OS*1E3,rot90(plot_varINTERP(:,1:end),3));
axis image

xlabel('X [mm]','Interpreter','latex')
ylabel('Z [mm]','Interpreter','latex' )
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
start_z = round(length(an.domain.zpar_OS)*0.5);
end_z = round(length(an.domain.zpar_OS)*0.7);
imagesc(an.domain.xpar_OS*1E3,an.domain.zpar_OS(start_z:end_z)*1E3,rot90(plot_varINTERP(:,start_z:end_z),3));
axis image;

xlabel('X [mm]','Interpreter','latex' )
ylabel('Z [mm]','Interpreter','latex' )
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
%%
end
%%
function plot_time_frame(an,frame,element_x)
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
p=imagesc(ax,an.domain.zpar*1E3,an.domain.xpar*1E3,squeeze(plot_value(1,:,:)) );
hold on;
caxis(ax,[max(max(max(plot_value)))-dBVALUES max(max(max(plot_value)))])
title_txt = 'Left-Aperture Pressure Field, $t_i$ = ';
title(ax,[title_txt,num2str(1)],'Interpreter','latex')
xlabel(ax,an.dslice.xlabel,'Interpreter','latex')
ylabel(ax,an.dslice.ylabel,'Interpreter','latex')

xlim(ax,[min(an.domain.zpar) max(an.domain.zpar)]*1E3)
ylim(ax,[min(an.domain.xpar) max(an.domain.xpar)]*1E3)

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
    set(p,'XData',an.domain.zpar*1E3,'YData',an.domain.xpar*1E3,'CData',squeeze(plot_value(i,:,:)))
    drawnow
end
plot([an.domain.zpar(1)*1E3 an.domain.zpar(end)*1E3],[0 0],'--','Color','white','LineWidth',2)
plot([0*1E3 an.domain.zpar(an.depth)*1E3],[element_x 0]*1E3,'--','Color','white','LineWidth',2)
%Capture the indexes of the max values in this plot
[max_val,max_idx] = max(squeeze(plot_value(i,:,:))');

% Plot the values in a figure
f3 = figure('WindowState','maximized');
plot(an.domain.zpar(max_idx)*1E3,an.domain.xpar*1E3,'LineWidth',2)

xlabel(an.dslice.xlabel,'Interpreter','latex')
ylabel(an.dslice.ylabel,'Interpreter','latex')
axis image;
xlim([an.domain.zpar(1) an.domain.zpar(end)]*1E3)

% Add 2 points in the same figure
hold on;
P_ZX = [an.domain.zpar(max_idx) ; an.domain.xpar]' * 1E3;
P1 = P_ZX(160,:);
P2 = P_ZX(191,:);
plot(P1(:,1),P1(:,2),'ro','LineWidth',2','MarkerSize',5)
plot(P2(:,1),P2(:,2),'ro','LineWidth',2','MarkerSize',5)

%Find many differences
P1_many = P_ZX(160,:);
P2_many = P_ZX(191,:);
Angle_many=90-atand((P2_many(:,2)-P1_many(:,2))./(P2_many(:,1)-P1_many(:,1)));
title(['Max pressure location in X-Z plane, ' num2str(an.theta) ' deg, @ Y = 0 mm '],'Interpreter','latex')
legend(['Angle between 2 points: ',num2str(mean(Angle_many)) ,' deg'],'Interpreter','latex')
plot([an.domain.zpar(1) an.domain.zpar(end)]*1E3,[0 0],'--','Color','black','LineWidth',2)

ax =gca;
ax.TickLabelInterpreter = 'latex';
set(findall(gca,'-property','FontSize'),'FontSize',40)

set(gcf, 'Color', 'w');
set(gca,'FontSize',18)
end
%%
function [an_press,an] = Compare_Signals(an, var, plot_val)
%%
[x_el,arr,td_x,apod_x,apod_right,apod_left] = xAM_Generate_TD_APOD(an.theta);
arr_length = arr.W/2*1E-3 ;
element_x = arr_length - an.domain.zpar(an.depth)/tand(90-an.theta);
td =  abs(sind(an.theta) * (arr_length - element_x)/an.medium.c0);

td_ofzstart = - (an.domain.start(3)-1)*an.domain.dxpar;
travel_dist_xwave = an.domain.zpar(an.depth)*cosd(an.theta) + (arr_length+an.domain.xpar(an.lat))*sind(an.theta); % m This is true only for the middle
an_td = travel_dist_xwave/an.medium.c0 ; % There is some small difference because the dist that cross x is not always a gridpoint in the computational domain

travel_dist_xwave = an.domain.zpar(an.depth) + arr_length*sind(an.theta); % m This is true only for the middle
c_x = an.medium.c0 / cosd(an.theta);
an_td = travel_dist_xwave/c_x ; % There is some small difference because the dist that cross x is not always a gridpoint in the computational domain

frame = round((an_td + an.gauss_dl)/an.domain.dtpar );

tdshift = an_td - an.td;
an_press = var;
an_press  = real(ifft(fft(an_press).*exp(-1i*an.omega*tdshift)));%;

if (plot_val)
    
    figure('WindowState','maximized');
    hold on;plot(an.domain.tpar*1E6*an.OSFactor,resampleSINC(squeeze(an.plin.data(:,an.lat,an.depth))',an.OSFactor),'LineWidth',2)
    hold on;plot(an.domain.tpar*1E6*an.OSFactor,an_press,'-.','LineWidth',2)
    legend('INCS','Analytical')
    
    ax = gca;
    xlabel('Time [$\mu$s]','Interpreter','latex')
    ylabel('Pressure [kPa]','Interpreter','latex')
    title(['Axial profile @ Z =', num2str(an.domain.zpar(an.depth)*1E3),' [mm]'],'Interpreter','latex')
    %     xlim( (tdshift - [0 1]*an.gauss_dl/2 + [0 1] * (an.gauss_win * 2 * an.medium.Fnyq)) * 1E6)
    
    ax.TickLabelInterpreter = 'latex';
    set(findall(gca,'-property','FontSize'),'FontSize',40)
    set(gcf,'Color','white')
    grid on;
    grid minor;
else
    close gcf
end
end