function xAM_Compare_Input_Signals
T_width = 1.5/medium.freq0;
T_delay = 3.0/medium.freq0 ;% /f0
Power = 2;
OSFactor = 1;
load('15MHz_pulse_for_Agis.mat')
domain.dtpar_new = 1/(2*medium.freq0*domain.Fnyq)/OSFactor;
domain.tpar_new = (1:(domain.tdimpar+1)*OSFactor)* domain.dtpar_new;
INCS_pulse=1E5*exp( - ((domain.tpar_new-T_delay)/(T_width/2)).^Power).* sin(2*pi*(domain.tpar_new*medium.freq0-T_delay)).* (1 + sign(domain.tpar_new))/2	;
% INCS_pulse =max(HF_pulse) /  max(INCS_pulse) * INCS_pulse;
max(INCS_pulse)

figure('WindowState','Maximized');
h_f2 = subplot(2,1,1);
ax = gca;   hold on;
plot(domain.tpar_new,INCS_pulse*1E-3,'LineWidth',2)
% plot([1:length(HF_pulse)]*1/sampling_freq + T_delay/2+T_width/2,HF_pulse*1E-3,'--','LineWidth',2)
% plot(t*1E-6,pulse*1E-3,'--','LineWidth',2)
% xlim([ 1/sampling_freq+ T_delay-T_width length(HF_pulse)*1/sampling_freq+ T_delay-T_width])
legend(['INCS - Tw = ' num2str(T_width*medium.freq0) ', Power = ' , num2str(Power)],'kWave')
ylabel('\textit{p} [kPa]')
xlabel('\textit{t} [$\mu$s]')
title('Time Signature'); box on;grid on;
set(gca,'FontSize',20,'XTickLabel',ax.XTick*1e6)

%%
% sampling_freq=1/(t(2)-t(1))*1E6;
F_s = 1/(domain.dtpar_new*2);
[f_kwave,Sp_psim] = Freq_Calc(HF_pulse,sampling_freq/2);            % Spectrum of simulation
[f_INCS,Sp_pan] = Freq_Calc(INCS_pulse,F_s);                % Spectrum of analytical pulse

% Create plot of frequency spectrum for driving pulse and radial oscillations
%                     f7 = figure('WindowState','maximized');
h_f2 = subplot(2,1,2);
ax = gca;   hold on;
m_plot=plot(f_INCS*1E-6,Sp_pan-max(Sp_pan),'-','LineWidth', 1, 'color', 'blue','LineWidth',2);
pr_plot=plot(f_kwave*1E-6,Sp_psim-max(Sp_psim), '--', 'LineWidth', 1,'color', 'red','LineWidth',2);
legend([m_plot , pr_plot], 'INCS', 'kWave')
ylabel('normalized amplitude [dB]');
xlabel('Frequency [MHz]');
title('Frequency spectrum');
box on; grid on;

ylim([-45 max(Sp_pan(1:domain.tdimpar)-max(Sp_pan))]);
xlim([0 F_s*1E-6/OSFactor]);

set(gca,'FontSize',20)
hold off;

%% Compare with INCS Results
T_width = 1.87/medium.freq0;
T_delay = 6.0/medium.freq0 ;% /f0
Power = 2;
INCS_pulse=resampleSINC(8E4*exp( - ((domain.tpar-T_delay)/(T_width/2)).^Power).* sin(2*pi*medium.freq0*(domain.tpar-T_delay)),OSFactor);
HF_pulse=resampleSINC(pnl.xwave_data(1:length_t,round(dims(1)/2)+1,depth_pos)',OSFactor);

figure('WindowState','Maximized');
h_f2 = subplot(2,1,1);
ax = gca;   hold on;
plot(domain.dtpar/OSFactor+domain.dtpar/OSFactor*[1:domain.tdimpar*OSFactor],INCS_pulse*1E-3,'LineWidth',2)
plot(domain.dtpar/OSFactor*[1:length_t*OSFactor],HF_pulse*1E-3,'--','LineWidth',2)
xlim([ 0 domain.dtpar*100])
legend(['Input - Tw = ' num2str(T_width*medium.freq0) ', Power = ' , num2str(Power)],['INCS @ ',num2str(round(domain.par{3}(depth_pos),2)),' [mm]'])
ylabel('Pressure [kPa]')
xlabel('Time [\mu sec]')
title('Time Signature'); box on;grid on;
set(gca,'FontSize',20,'XTickLabel',ax.XTick*1e6)

%%
F_s = 1/(domain.dtpar_new*2);
[f_INCS,Sp_psim] = Freq_Calc(HF_pulse,F_s);                  % Spectrum of simulation
[f_Input,Sp_pan] = Freq_Calc(INCS_pulse,F_s);                % Spectrum of analytical pulse

% Create plot of frequency spectrum for driving pulse and radial oscillations
%                     f7 = figure('WindowState','maximized');
h_f2 = subplot(2,1,2);
ax = gca;   hold on;
m_plot=plot(f_Input*1E-6,Sp_pan-max(Sp_pan),'-','LineWidth', 1, 'color', 'blue','LineWidth',2);
pr_plot=plot(f_INCS*1E-6,Sp_psim-max(Sp_psim), '--', 'LineWidth', 1,'color', 'red','LineWidth',2);
legend([m_plot , pr_plot], ['Input - Tw = ' num2str(T_width*medium.freq0) ', Power = ' , num2str(Power)],['INCS @ ',num2str(round(domain.par{3}(depth_pos),2)),' [mm]'])
ylabel('normalized amplitude [dB]');
xlabel('Frequency [MHz]');
title('Frequency spectrum');
box on; grid on;

ylim([-80 max(Sp_pan(1:domain.tdimpar)-max(Sp_pan))]);
xlim([0 F_s*1E-6/OSFactor]);

set(gca,'FontSize',20)
hold off;
end

% Calculate the frequency spectrum based on a sampling frequency and the pressure values
function [f_p,S_p] = Freq_Calc(Pres,F_s)

Nfft_p = 2^nextpow2(2*length(Pres));                      % Zero-padding to improve fft performance until length equal to 2N
Y_p = abs(fft(Pres,Nfft_p));                              % Convert the driving pulse to the frequency domain
f_p = linspace(0,1,Nfft_p/2+1)*F_s;                     % [Hz] Frequency sample Number
% Define the frequency domain for half of sampling Frequency [0,Fs/2] - Nyquist, equally spaced values between  0 and 1
% Other way to define frequency domain : f_p = Fs*(0:(Nfft_p/2))/Nfft_p
Np = 1:Nfft_p/2+1;                                         % Creating a length vector in order to get the one-sided fft values , It is equal with [1:length(f_p)]
S_p = 20*log10(Y_p(Np));                                     % Convert pressure to amplitude [db]

end