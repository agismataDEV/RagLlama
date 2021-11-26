syms x y z t w c0 rho P0 theta beta_mn
clc
% p_list = [ sin(w*t - w/c0*sind(theta)*x-w/c0*cosd(theta)*z) ; sin(x+c0*t) ; cos(x-c0*t) ; cos(x+c0*t); sin(2*(x+c0*t)); cos(x+c0*t) - sin(x+c0*t) ; sin(x-c0*t) + cos(x-c0*t) ; ...
%     cos(x+c0*t) + cos(x-c0*t); sin(x+c0*t) - cos(x-c0*t); P0* (sin(w*t - w/c0*sind(theta)*x - w/c0*cosd(theta)*z) - sin(w*t + w/c0*sind(theta)*x - w/c0*cosd(theta)*z))];
p_list = [P0* (sin(w*t - w/c0*sind(theta)*x - w/c0*cosd(theta)*z) + sin(w*t + w/c0*sind(theta)*x - w/c0*cosd(theta)*z))];
for i =1:length(p_list)
    disp(' ')
    p_lagr = p_list(i);
%     p_lagr = vpa(simplify(subs(p_lagr, [rho c0 P0 w theta], [1060 1480 4E4 2*pi*freq0 21])));
    p_string  = strjoin({char("By using an input pressure equal to : ") , char(p_lagr)});
    disp(p_string);
    phi_lagr = -1/rho * int(p_lagr,t);
%     phi_sq = simplify(phi_lagr^2);
%     phi_sq_x = simplify(diff(diff(phi_sq,x),x));
%     phi_sq_y = simplify(diff(diff(phi_sq,y),y));
%     phi_sq_z = simplify(diff(diff(phi_sq,z),z));
%     phi_sq_t = simplify(diff(diff(phi_sq,t),t));
%     Lagrangian_simpl = subs(rho*simplify(phi_sq_x + phi_sq_y + phi_sq_z - 1/c0^2 * phi_sq_t)/4 , theta, 0  );
%     L_String = strjoin({char("The pressure due to Lagrangian equals to : "),char(Lagrangian_simpl)});
%     disp(L_String)
    
    phi_sq_x = simplify(diff(phi_lagr,x)^2);
    phi_sq_y = simplify(diff(phi_lagr,y)^2);
    phi_sq_z = simplify(diff(phi_lagr,z)^2);
    potential = p_lagr^2/(2*rho*c0^2);
    Lagrangian = rho*simplify(phi_sq_x + phi_sq_y + phi_sq_z)/2 - potential;
%     Lagrangian_sub = subs(Lagrangian , theta, 21  );
%     Lagrangian_sub = subs(Lagrangian_sub , [rho c0 P0 w theta x z ], [1060 1480 4E4 2*pi*freq0 21 0 domain.parz(143)*1E-3]  );
%     L_String = strjoin({char("The Lagrangian equals to : "),char(Lagrangian_sub)});
%     disp(L_String)
    
    L1_X = diff(diff(Lagrangian,x),x);
    L1_Y = diff(diff(Lagrangian,y),y);
    L1_Z = diff(diff(Lagrangian,z),z);
    L2_T = diff(diff(Lagrangian,t),t)/c0^2;
    Lagrangian_NL = simplify(L1_X+L1_Y+L1_Z+L2_T);
    L_String = strjoin({char("The nonlinear term due to Lagrangian equals to : "),char(Lagrangian_NL)});
%     Lagrangian_NL_sub = vpa(simplify(subs(Lagrangian_NL , [P0 theta c0 rho x y z w ],[1E6 15 1480 1060 0 0 domain.par{3}(300)*1E-3 2*pi*8*15E6] )));
%     Lagrangian_NL_sub_zt = subs(Lagrangian_NL_sub, t ,domain.tpar);
%     L_String = strjoin({char("The nonlinear term due to Lagrangian equals to : "),char(Lagrangian_NL_sub)});
%     figure;plot(domain.part+domain.parz(142)*1E-3/1480+max(td_x)*1E-6,real(double(subs(Lagrangian_sub * phi.blackman_win * phi.sign_win,t,domain.part))))
%     disp(L_String)
    
    Medium_NL = diff(diff(p_lagr^2,t),t)/c0^4/rho*beta_mn;
    L_String = strjoin({char("The nonlinear term due to Medium nonlinearities equals to : "),char(Medium_NL)});
%     Medium_NL_sub = vpa(simplify(subs(Medium_NL , [P0 theta c0 rho x y z w beta_mn ],[1E6 15 1480 1060 0 0 domain.par{3}(300)*1E-3 2*pi*8*15E6 3] )));
%     Medium_NL_sub_zt = subs(Medium_NL_sub, t ,domain.tpar);
%     L_String = strjoin({char("The nonlinear term due to Medium nonlinearities equals to : "),char(Medium_NL_sub)});
    disp(L_String)
end

%% Computed from matlab symbolic expressions
freq0 = 15E6;
t_interest = domain.tpar(50:100);
p_lag_an_expr = vpa(simplify(subs(Lagrangian_NL, [rho c0 P0 w theta y ], [ 1060 1480 5E4 2*pi*freq0 21 0 ])));
p_lag_an_expr = subs(p_lag_an_expr,t, t_interest);
%% Generate the result for all x and t , (t,x)
x_interest = domain.par{1}(75:140);
p_lag_an_x = reshape(subs(p_lag_an_expr, x, x_interest*1E-3)', [length(x_intetest) length(t_interst)])';
%% Generate the result for (t,x,z)
z_interest = domain.par{3}(90:180);
p_lag_an_xz=subs(p_lag_an_x, z, z_interest*1E-3);
%% Reshape the results in order to match the simulation results
p_lag_an_xz_expr=real(reshape(p_lag_an_xz',[length(z_intetest),length(x_intetest),length(t_interst)]));
p_lag_an_xz=permute(p_lag_an_xz_expr,[3 2 1]);
p_tot = double(real(p_lag_an_xz));

%% Computed from analytical expressions on paper  
rho = 1060;
c0 = 1480;
P0 = 5E4;
freq0 = 15E6;
w = 2* pi * freq0;
theta = 21;
t = domain.tpar;
x = domain.par{1}(75:140);
z = domain.par{3}(90:180);
[X,T,Z]=meshgrid(x,t,z);
blackman_L             = 3;
blackman_delay         = 35.5; 
blwin                  = 0.35869 - 0.48829 * cos(2*pi*(freq0*T-blackman_delay)/blackman_L) + 0.14128*cos(4*pi*(freq0*T-blackman_delay)/blackman_L) - 0.01168*sin(6*pi*(freq0*T-blackman_delay)/blackman_L);
signature              = (sign(T-blackman_delay/freq0)-sign(T-blackman_delay/freq0-blackman_L/freq0))/2;
window                 = 1;blwin .* signature;
% window(:,1:20,:)            = 0;
% window(:,40:end,:)            = 0; 
% window = zeros(size(X));
% window(30:40,:,:) = 1;

N2 = 8*w^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*T - 2* w/c0 * cosd(theta) * Z).*window-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*T - 2 * w/c0* cosd(theta) * Z).*window + cos(2*w/c0*sind(theta) * X).*window);
figure('WindowState','Maximized');imagesc(z,x,squeeze(max(abs(N2/freq0^2*c0^2))));colorbar;
xlabel('Depth [mm]')
ylabel('X [mm]')
title('Analytical for Source Term due to Lagrangian, @ 21 deg')

xAM_Colormaps(file)
set(gcf,'Color','white')
set(gca,'FontSize',20)

