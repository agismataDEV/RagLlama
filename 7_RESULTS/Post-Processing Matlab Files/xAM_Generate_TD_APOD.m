%% Phased Array Specs
% clear all;
% close all;
function [x_el,arr,td_x,apod_x,apod_right,apod_left] = xAM_Generate_TD_APOD(theta_list)
for i = 1:length(theta_list)
    i_theta = theta_list(i);
location = '';
% mkdir('../../5_Input_Files/', location)
c0 = 1480;
freq0 = 15E6;

theta = i_theta; % deg
arr.N_el = 64;
arr.H_el = 4;    %mm
arr.W_el = 0.1; %mm
arr.Kerf = 0.00; %mm

arr.W = (arr.W_el + arr.Kerf)* (arr.N_el-1) + arr.W_el;
arr.H = arr.H_el ;

x_el = arr.W_el/2 + (0:arr.N_el-1)*( arr.W_el + arr.Kerf);
%% XWAVE
td_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_x.dat'];
td_x = abs(x_el - arr.W/2 )*sind(i_theta);
td_x = abs(td_x - max(td_x))'*1E-3/c0*1E6; %usec
size_td = size(td_x);
save(td_filename, 'size_td','td_x','-ASCII', '-DOUBLE');

apod_x_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_x.dat'];
% If even , all ones, if odd, the middle is silenced
% apod_x = [ones(floor(arr.N_el/2),1); zeros(arr.N_el - 2*floor(arr.N_el/2),1) ; ones(floor(arr.N_el/2),1)];
apod_x = [tukeywin(arr.N_el/2,0.5); tukeywin(arr.N_el/2,0.5)];
if (i_theta==0) ;apod_x = tukeywin(arr.N_el,0.2);end

size_apod = size(apod_x);
save(apod_x_filename, 'size_apod','apod_x','-ASCII', '-DOUBLE');
%% LEFT
td_left_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_left.dat'];
td_left = abs(x_el)*sind(i_theta);
td_left =abs(td_left - min(td_left))'*1E-3/c0*1E6; %usec
size_td = size(td_left);
save(td_left_filename, 'size_td','td_left','-ASCII', '-DOUBLE');

apod_left_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_left.dat'];
apod_left = [ones(floor(arr.N_el/2),1);zeros(arr.N_el - floor(arr.N_el/2),1)];
apod_left = apod_left.*apod_x;
size_apod = size(apod_left);
save(apod_left_filename, 'size_apod','apod_left','-ASCII', '-DOUBLE');
%% RIGHT
td_right_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_right.dat'];
td_right = abs(x_el)*sind(i_theta);
td_right = abs(td_right - max(td_right))'*1E-3/c0*1E6; %usec
size_td = size(td_right);   
save(td_right_filename, 'size_td','td_right','-ASCII', '-DOUBLE');

apod_right_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_right.dat'];
apod_right = [zeros(arr.N_el - floor(arr.N_el/2),1);ones(floor(arr.N_el/2),1)];
apod_right = apod_right.*apod_x;
size_apod = size(apod_right);
save(apod_right_filename, 'size_apod','apod_right','-ASCII', '-DOUBLE');
end
end