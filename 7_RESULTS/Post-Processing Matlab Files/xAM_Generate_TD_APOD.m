%% Phased Array Specs
% clear all;
% close all;
function [arr] = xAM_Generate_TD_APOD(theta_list)
for i = 1:length(theta_list)
    i_theta = theta_list(i);
location = '';
% mkdir('../../5_Input_Files/', location)
c0 = 1480;

arr.N_el = 22;
arr.W_el = 0.08; %mm
arr.Kerf = 0.02; %mm

arr.W = (arr.W_el + arr.Kerf)* (arr.N_el-1) + arr.W_el;

arr.x_el = arr.W_el/2 + (0:arr.N_el-1)*( arr.W_el + arr.Kerf);
arr.apod_val = 0.3;
%% XWAVE
td_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_x.dat'];
arr.td_x = abs(arr.x_el - arr.W/2 )*sind(i_theta);
arr.td_x = abs(arr.td_x - max(arr.td_x))'*1E-3/c0*1E6; %usec
td_x=arr.td_x;
size_td = size(arr.td_x);
save(td_filename, 'size_td','td_x','-ASCII', '-DOUBLE');

apod_x_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_x.dat'];
% If even , all ones, if odd, the middle is silenced
% apod_x = [ones(floor(arr.N_el/2),1); zeros(arr.N_el - 2*floor(arr.N_el/2),1) ; ones(floor(arr.N_el/2),1)];
arr.apod_x = [tukeywin(arr.N_el/2,arr.apod_val); tukeywin(arr.N_el/2,arr.apod_val)];
if (i_theta==0) ;arr.apod_x = tukeywin(arr.N_el,arr.apod_val);end

size_apod = size(arr.apod_x);
apod_x=arr.apod_x;
save(apod_x_filename, 'size_apod','apod_x','-ASCII', '-DOUBLE');
%% LEFT
td_left_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_left.dat'];
arr.td_left = abs(arr.x_el)*sind(i_theta);
arr.td_left =abs(arr.td_left - min(arr.td_left))'*1E-3/c0*1E6; %usec
size_td = size(arr.td_left);
td_left=arr.td_left;
save(td_left_filename, 'size_td','td_left','-ASCII', '-DOUBLE');

apod_left_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_left.dat'];
arr.apod_left = [ones(floor(arr.N_el/2),1);zeros(arr.N_el - floor(arr.N_el/2),1)];
arr.apod_left = arr.apod_left.*apod_x;
if (i_theta==0) ;arr.apod_left =apod_x;end
size_apod = size(arr.apod_left);
apod_left=arr.apod_left;
save(apod_left_filename, 'size_apod','apod_left','-ASCII', '-DOUBLE');
%% RIGHT
td_right_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_right.dat'];
arr.td_right = abs(arr.x_el)*sind(i_theta);
arr.td_right = abs(arr.td_right - max(arr.td_right))'*1E-3/c0*1E6; %usec
size_td = size(arr.td_right);   
td_right=arr.td_right;
save(td_right_filename, 'size_td','td_right','-ASCII', '-DOUBLE');

apod_right_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_right.dat'];
arr.apod_right = [zeros(arr.N_el - floor(arr.N_el/2),1);ones(floor(arr.N_el/2),1)];
arr.apod_right = arr.apod_right.*arr.apod_x;
if (i_theta==0) ;arr.apod_right =arr.apod_x;end
size_apod = size(arr.apod_right);
apod_right=arr.apod_right;
save(apod_right_filename, 'size_apod','apod_right','-ASCII', '-DOUBLE');


%% PAM
idx_s=(64-arr.N_el)/2+1; idx_e=(64-arr.N_el)/2+arr.N_el;
td_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_pam.dat'];
arr.td_pam = zeros(64,1);
% arr.td_pam(idx_s:idx_e) = sqrt((arr.x_el - arr.W/2).^2+ (4.93)^2);
% arr.td_pam(idx_s:idx_e) = abs(arr.td_pam(idx_s:idx_e) - max(arr.td_pam(idx_s:idx_e)))'*1E-3/c0*1E6; %usec
arr.td_pam(idx_s:idx_e) = sqrt((arr.x_el - arr.W/2).^2+ (3)^2);
arr.td_pam(idx_s:idx_e) = abs(arr.td_pam(idx_s:idx_e) - min(arr.td_pam(idx_s:idx_e)))'*1E-3/c0*1E6; %usec
% arr.td_pam = arr.td_x;
size_td = size(arr.td_pam);
td_pam=arr.td_pam;
save(td_filename, 'size_td','td_pam','-ASCII', '-DOUBLE');

apod_pam_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_pam.dat'];
% If even , all ones, if odd, the middle is silenced
% apod_x = [ones(floor(arr.N_el/2),1); zeros(arr.N_el - 2*floor(arr.N_el/2),1) ; ones(floor(arr.N_el/2),1)];
arr.apod_pam = zeros(64,1);
arr.apod_pam(idx_s:idx_e) = [tukeywin(arr.N_el/2,arr.apod_val); tukeywin(arr.N_el/2,arr.apod_val)];
arr.apod_pam(idx_s + arr.N_el/4:idx_e - arr.N_el/4)=1;
% arr.apod_pam=arr.apod_x;
if (i_theta==0) ;arr.apod_pam = tukeywin(arr.N_el,arr.apod_val);end

size_apod = size(arr.apod_pam);
apod_pam=arr.apod_pam;
save(apod_pam_filename, 'size_apod','apod_pam','-ASCII', '-DOUBLE');
%% ODD
td_odd_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_odd.dat'];
save(td_odd_filename, 'size_td','td_pam','-ASCII', '-DOUBLE');

apod_odd_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_odd.dat'];
arr.apod_odd = ones(64,1);
arr.apod_odd(2:2:64) = 0;
arr.apod_odd = arr.apod_odd.*arr.apod_pam;
size_apod = size(arr.apod_odd);
apod_odd=arr.apod_odd;
save(apod_odd_filename, 'size_apod','apod_odd','-ASCII', '-DOUBLE');
%% EVEN
td_even_filename = ['../../5_Input_Files/TD/', location,'td_' num2str(i_theta*10) 'deg_even.dat'];
save(td_even_filename, 'size_td','td_pam','-ASCII', '-DOUBLE');

apod_even_filename = ['../../5_Input_Files/APOD/', location,'apod_' num2str(i_theta*10) 'deg_even.dat'];
arr.apod_even = ones(64,1);
arr.apod_even(1:2:64) = 0;
arr.apod_even = arr.apod_even.*arr.apod_pam;
size_apod = size(arr.apod_even);
apod_even=arr.apod_even;
save(apod_even_filename, 'size_apod','apod_even','-ASCII', '-DOUBLE');
end
end