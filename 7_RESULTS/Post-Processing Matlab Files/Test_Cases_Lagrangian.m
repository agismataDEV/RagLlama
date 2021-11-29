syms x z y z t w c0 rho P0 theta beta_mn
clc
% p_list = [ sin(w*t - w/c0*sind(theta)*x-w/c0*cosd(theta)*z) ; sin(x+c0*t) ; cos(x-c0*t) ; cos(x+c0*t); sin(2*(x+c0*t)); cos(x+c0*t) - sin(x+c0*t) ; sin(x-c0*t) + cos(x-c0*t) ; ...
%     cos(x+c0*t) + cos(x-c0*t); sin(x+c0*t) - cos(x-c0*t); P0* (sin(w*t - w/c0*sind(theta)*x - w/c0*cosd(theta)*z) - sin(w*t + w/c0*sind(theta)*x - w/c0*cosd(theta)*z))];
blackman_L             = 6;
blackman_delay         = 6; 
blwin                  = 0.35869 - 0.48829 * cos(2*pi*(freq0*t-blackman_delay)/blackman_L) + 0.14128*cos(4*pi*(freq0*t-blackman_delay)/blackman_L) - 0.01168*sin(6*pi*(freq0*t-blackman_delay)/blackman_L);
signature              = (sign(t-blackman_delay/freq0)-sign(t-blackman_delay/freq0-blackman_L/freq0))/2;

p_list = [P0* (sin(w*t - w/c0*sind(theta)*x - w/c0*cosd(theta)*z) + sin(w*t + w/c0*sind(theta)*x - w/c0*cosd(theta)*z))*blwin*signature];
for i =1:length(p_list)
    disp(' ')
    p_lagr = p_list(i);
%     p_lagr = vpa(simplify(subs(p_lagr, [rho c0 P0 w theta], [1060 1480 4E4 2*pi*freq0 21])));
    p_string  = strjoin({char("By using an input pressure equal to : ") , char(p_lagr)});
    disp(p_string);
%     phi_lagr = -1/rho * int(p_lagr,t);
% %     phi_sq = simplify(phi_lagr^2);
% %     phi_sq_xyz = diff(phi_sq,x,2) +  diff(phi_sq,y,2) + diff(phi_sq,z,2);
% %     phi_sq_t = diff(phi_sq,t,2)/c0^2;
% %     Lagrangian_simpl = simplify(phi_sq_xyz - phi_sq_t);
%     
%     grad_phi    = simplify(diff(phi_lagr,x,2)+ diff(phi_lagr,y,2) + diff(phi_lagr,z,2));
%     potential   = p_lagr^2/(2*rho*c0^2);
%     Lagrangian  = rho*grad_phi/2 - potential;
%     
%     L2_XYZ  =   diff(Lagrangian,x,2) + diff(Lagrangian,y,2)+diff(Lagrangian,z,2);
%     L2_T    =   diff(Lagrangian,t,2)/c0^2;
%     Lagrangian_NL = simplify(L2_XYZ+L2_T);
%     L_String = strjoin({char("The nonlinear term due to Lagrangian equals to : "),char(Lagrangian_NL)});
%     Lagrangian_NL_sub = vpa(simplify(subs(Lagrangian_NL , [P0 theta c0 rho x y z w ],[1E6 15 1480 1060 0 0 domain.par{3}(300)*1E-3 2*pi*8*15E6] )));
%     Lagrangian_NL_sub_zt = subs(Lagrangian_NL_sub, t ,domain.tpar);
%     L_String = strjoin({char("The nonlinear term due to Lagrangian equals to : "),char(Lagrangian_NL_sub)});
%     figure;plot(domain.part+domain.parz(142)*1E-3/1480+max(td_x)*1E-6,real(double(subs(Lagrangian_sub * phi.blackman_win * phi.sign_win,t,domain.part))))
%     disp(L_String)
    
    Medium_NL = diff(p_lagr^2,t,2)/c0^4/rho*beta_mn;
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
p_str =  char(p_lag_an_expr);
p_str = strrep(char(p_str),'*','.*');
p_str = strrep(char(p_str),'^','.^');
p_str = strrep(char(p_str),'/','./')
%%
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

%% Computed from analytical expressions on paper  , Lagrangian
rho = 1060;
c0 = 1480;
P0 = 5E4;
freq0 = 1/(an.domain.tpar(2)-an.domain.tpar(1));
w = 2* pi * freq0;
theta = 21;

T = an.domain.tpar(30:60);
X = an.domain.xpar*1E3;
Z = an.domain.zpar*1E3;
[x,t,z]=meshgrid(X,T,Z);
blackman_L             = 6;
blackman_delay         = 6; 
blwin                  = 0.35869 - 0.48829 * cos(2*pi*(freq0*t-blackman_delay)/blackman_L) + 0.14128*cos(4*pi*(freq0*t-blackman_delay)/blackman_L) - 0.01168*sin(6*pi*(freq0*t-blackman_delay)/blackman_L);
signature              = (sign(t-blackman_delay/freq0)-sign(t-blackman_delay/freq0-blackman_L/freq0))/2;
window                 = blwin ;%.* signature;
% window(:,1:20,:)            = 0;
% window(:,40:end,:)            = 0; 
% window = zeros(size(X));
% window(30:40,:,:) = 1;

% N2 = 8*w.^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*t - 2 * w/c0* cosd(theta) * z) + cos(2*w/c0*sind(theta) * x));
% N2 = 4*w^2*P0^2/rho/c0^4*(sind(theta)^2*(1+cosd(theta)^2)*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-sind(theta)^2*cos(2*w/c0*sind(theta) * x));
N2 = 8732921017.4931237281138599603576.*cos(45642.410171900795453546919281243.*x) - 8732921017.4931237281138599603576.*cos(188495559.2153875827789306640625.*t - 118902.54362508789662473735720873.*z) + 17448838406227485.574047986159216.*cos(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).*cos(22821.205085950397726773459640622.*x) - 15222746082.945919248336281588498.*cos(22821.205085950397726773459640622.*x - 94247779.60769379138946533203125.*t + 59451.271812543948312368678604365.*z).*cos(94247779.60769379138946533203125.*t + 22821.205085950397726773459640622.*x - 59451.271812543948312368678604365.*z) - 17465842034.986247456227719920715.*cos(22821.205085950397726773459640622.*x - 94247779.60769379138946533203125.*t + 59451.271812543948312368678604365.*z).^2 - 17465842034.986247456227719920715.*cos(94247779.60769379138946533203125.*t + 22821.205085950397726773459640622.*x - 59451.271812543948312368678604365.*z).^2 + 17465842034.986247456227719920715;

N2 = N2/freq0^2*c0^2;
figure('WindowState','Maximized');imagesc(Z,X,squeeze(max(abs(N2))));colorbar;
xlabel('Depth [mm]')
ylabel('X [mm]')
title(['Analytical for Source Term due to Lagrangian, @ ', num2str(theta), ' deg'])

xAM_Colormaps(file)
set(gcf,'Color','white')
% set(gca,'FontSize',20)
xlim([1.94 2.15])
ylim([-0.2 0.2])

figure('WindowState','Maximized');plot(Z,squeeze(max(abs(N2(:,length(X)/2+1,:)))))
xlabel('Depth [mm]')
ylabel('Pressure')
grid on;
set(gcf,'Color','white')
% set(gca,'FontSize',20)
%% Computed from analytical expressions on paper  
rho = 1060;
c0 = 1480;
P0 = 5E4;
freq0 = 1/(an.domain.tpar(2)-an.domain.tpar(1));
w = 2* pi * freq0;
theta = 21;

T = an.domain.tpar(30:60);
X = an.domain.xpar*1E3;
Z = an.domain.zpar*1E3;
[x,t,z]=meshgrid(X,T,Z);
blackman_L             = 6;
blackman_delay         = 6; 
blwin                  = 0.35869 - 0.48829 * cos(2*pi*(freq0*t-blackman_delay)/blackman_L) + 0.14128*cos(4*pi*(freq0*t-blackman_delay)/blackman_L) - 0.01168*sin(6*pi*(freq0*t-blackman_delay)/blackman_L);
signature              = (sign(t-blackman_delay/freq0)-sign(t-blackman_delay/freq0-blackman_L/freq0))/2;
window                 = blwin ;%.* signature;
% window(:,1:20,:)            = 0;
% window(:,40:end,:)            = 0; 
% window = zeros(size(X));
% window(30:40,:,:) = 1;

% N2 = 8*w.^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*t - 2 * w/c0* cosd(theta) * z) + cos(2*w/c0*sind(theta) * x));
% N2 = 4*w^2*P0^2/rho/c0^4*(sind(theta)^2*(1+cosd(theta)^2)*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-sind(theta)^2*cos(2*w/c0*sind(theta) * x));
N2 = 8732921017.4931237281138599603576.*cos(45642.410171900795453546919281243.*x) - 8732921017.4931237281138599603576.*cos(188495559.2153875827789306640625.*t - 118902.54362508789662473735720873.*z) + 17448838406227485.574047986159216.*cos(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).*cos(22821.205085950397726773459640622.*x) - 15222746082.945919248336281588498.*cos(22821.205085950397726773459640622.*x - 94247779.60769379138946533203125.*t + 59451.271812543948312368678604365.*z).*cos(94247779.60769379138946533203125.*t + 22821.205085950397726773459640622.*x - 59451.271812543948312368678604365.*z) - 17465842034.986247456227719920715.*cos(22821.205085950397726773459640622.*x - 94247779.60769379138946533203125.*t + 59451.271812543948312368678604365.*z).^2 - 17465842034.986247456227719920715.*cos(94247779.60769379138946533203125.*t + 22821.205085950397726773459640622.*x - 59451.271812543948312368678604365.*z).^2 + 17465842034.986247456227719920715;

N2 = N2/freq0^2*c0^2;
figure('WindowState','Maximized');imagesc(Z,X,squeeze(max(abs(N2))));colorbar;
xlabel('Depth [mm]')
ylabel('X [mm]')
title(['Analytical for Source Term due to Lagrangian, @ ', num2str(theta), ' deg'])

xAM_Colormaps(file)
set(gcf,'Color','white')
% set(gca,'FontSize',20)
xlim([1.94 2.15])
ylim([-0.2 0.2])

figure('WindowState','Maximized');plot(Z,squeeze(max(abs(N2(:,length(X)/2+1,:)))))
xlabel('Depth [mm]')
ylabel('Pressure')
grid on;
set(gcf,'Color','white')
% set(gca,'FontSize',20)
%% Medium NL
%% Computed from matlab symbolic expressions
freq0 = 15E6;
t_interest = domain.tpar(50:100);
p_lag_an_expr = vpa(simplify(subs(Medium_NL, [rho c0 P0 w theta y beta_mn x ], [ 1060 1480 5E4 2*pi*freq0 21 0 3.21 0])));
p_str =  char(p_lag_an_expr);
p_str = strrep(char(p_str),'*','.*');
p_str = strrep(char(p_str),'^','.^');
p_str = strrep(char(p_str),'/','./');
%% Computed from analytical expressions on paper , MediumNL
rho = 1060;
c0 = 1480;
P0 = 5E4;
freq0 = 1/(an.domain.tpar(2)-an.domain.tpar(1));
w = 2* pi * freq0;
theta = 21;

T = an.domain.tpar;
Z = an.domain.zpar*1E3;
[z,t]=meshgrid(Z,T);

% N2 = 8*w.^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*t - 2 * w/c0* cosd(theta) * z) + cos(2*w/c0*sind(theta) * x));
MediumNL_an = 0.000012623573156503650292503715851399.*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).^2.*(4403256.2632714542030292409660046.*cos(376991118.43077518861551720599354.*t) - 61360331.072854405616356965504862.*sin(125663706.14359172953850573533118.*t) + 35507536.807933279098400180575178.*sin(251327412.28718345907701147066236.*t)).^2 + 0.000012623573156503650292503715851399.*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(1.0.*dirac(t - 0.000000049999999999999997737405591294313) - 1.0.*dirac(t - 0.000000099999999999999995474811182588626)).^2.*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869).^2 + 112130705864.61170866898196189099.*cos(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).^2.*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869).^2 - 112130705864.61170866898196189099.*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).^2.*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869).^2 - 0.000012623573156503650292503715851399.*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).^2.*(1659988503428021.4048190107102512.*sin(376991118.43077518861551720599354.*t) + 7710766612812676.6719851097765272.*cos(125663706.14359172953850573533118.*t) - 8924017342629789.34828119608616.*cos(251327412.28718345907701147066236.*t)).*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869) + 0.000012623573156503650292503715851399.*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).*(1.0.*dirac(1, t - 0.000000049999999999999997737405591294313) - 1.0.*dirac(1, t - 0.000000099999999999999995474811182588626)).*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869).^2 + 4758.974962863021911276636922599.*cos(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).*(1.0.*dirac(t - 0.000000049999999999999997737405591294313) - 1.0.*dirac(t - 0.000000099999999999999995474811182588626)).*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869).^2 + 0.000050494292626014601170014863405596.*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).^2.*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).*(1.0.*dirac(t - 0.000000049999999999999997737405591294313) - 1.0.*dirac(t - 0.000000099999999999999995474811182588626)).*(4403256.2632714542030292409660046.*cos(376991118.43077518861551720599354.*t) - 61360331.072854405616356965504862.*sin(125663706.14359172953850573533118.*t) + 35507536.807933279098400180575178.*sin(251327412.28718345907701147066236.*t)).*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869) + 4758.974962863021911276636922599.*cos(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).*sin(94247779.60769379138946533203125.*t - 59451.271812543948312368678604365.*z).*(0.5.*sign(t - 0.000000049999999999999997737405591294313) - 0.5.*sign(t - 0.000000099999999999999995474811182588626)).^2.*(4403256.2632714542030292409660046.*cos(376991118.43077518861551720599354.*t) - 61360331.072854405616356965504862.*sin(125663706.14359172953850573533118.*t) + 35507536.807933279098400180575178.*sin(251327412.28718345907701147066236.*t)).*(0.01168.*sin(376991118.43077518861551720599354.*t) + 0.48829.*cos(125663706.14359172953850573533118.*t) - 0.14128.*cos(251327412.28718345907701147066236.*t) - 0.35869);

MediumNL_an = MediumNL_an/freq0^2*c0^2;

figure('WindowState','Maximized');plot(Z,squeeze(max(abs(MediumNL_an))))
xlabel('Depth [mm]')
ylabel('Pressure')
grid on;
set(gcf,'Color','white')
% set(gca,'FontSize',20)

