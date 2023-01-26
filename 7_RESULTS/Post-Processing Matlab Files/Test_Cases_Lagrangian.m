syms x z y z t w c0 rho P0 theta beta_mn
clc
% p_list = [ sin(w*t - w/c0*sind(theta)*x-w/c0*cosd(theta)*z) ; sin(x+c0*t) ; cos(x-c0*t) ; cos(x+c0*t); sin(2*(x+c0*t)); cos(x+c0*t) - sin(x+c0*t) ; sin(x-c0*t) + cos(x-c0*t) ; ...
%     cos(x+c0*t) + cos(x-c0*t); sin(x+c0*t) - cos(x-c0*t); P0* (sin(w*t - w/c0*sind(theta)*x - w/c0*cosd(theta)*z) - sin(w*t + w/c0*sind(theta)*x - w/c0*cosd(theta)*z))];
% blackman_L             = 6;
% blackman_delay         = 6; 
% blwin                  = 0.35869 - 0.48829 * cos(2*pi*(freq0*t-blackman_delay)/blackman_L) + 0.14128*cos(4*pi*(freq0*t-blackman_delay)/blackman_L) - 0.01168*sin(6*pi*(freq0*t-blackman_delay)/blackman_L);
% signature              = (sign(t-blackman_delay/freq0)-sign(t-blackman_delay/freq0-blackman_L/freq0))/2;
% 
% p_list = [P0* (sin(w*t - w/c0*sind(theta)*x - w/c0*cosd(theta)*z) + sin(w*t + w/c0*sind(theta)*x - w/c0*cosd(theta)*z))*blwin*signature];

gauss_dl                  = 10.0/an.medium.freq0;
gauss_win                 = 5.0/an.medium.freq0;
% p_an  = exp( - ((t-gauss_dl - z*cosd(theta)/c0 - x*sind(theta)/c0)/(gauss_win/2)).^2).* sin(w*(t-gauss_dl) - w/c0*cosd(theta)*z - w/c0*sind(theta)*x);
p_an                =  P0*(sin(w*(t-gauss_dl - sind(theta)*x/c0 - cosd(theta)*z/c0))+ sin(w*(t-gauss_dl + sind(theta)*x/c0 - cosd(theta)*z/c0)));
p_list = [p_an];
%%
for i =1:length(p_list)
    disp(' ')
    p_lagr = p_list(i);
%     p_lagr = vpa(simplify(subs(p_lagr, [rho c0 P0 w theta], [1060 1480 4E4 2*pi*freq0 21])));
    p_string  = strjoin({char("By using an input pressure equal to : ") , char(p_lagr)});
    disp(p_string);
    phi_lagr = -1/rho * int(p_lagr,t);
% %     phi_sq = simplify(phi_lagr^2);
% %     phi_sq_xyz = diff(phi_sq,x,2) +  diff(phi_sq,y,2) + diff(phi_sq,z,2);
% %     phi_sq_t = diff(phi_sq,t,2)/c0^2;
% %     Lagrangian_simpl = simplify(phi_sq_xyz - phi_sq_t);
%     
    grad_phi    = simplify(diff(phi_lagr,x,2)+ diff(phi_lagr,y,2) + diff(phi_lagr,z,2));
    potential   = p_lagr^2/(2*rho*c0^2);
    Lagrangian  = simplify(rho*grad_phi/2 - potential);
%     
    L2_XYZ  =   diff(Lagrangian,x,2) + diff(Lagrangian,y,2)+diff(Lagrangian,z,2);
    L2_T    =   diff(Lagrangian,t,2)/c0^2;
    Lagrangian_NL = L2_XYZ+L2_T;
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
p_lag_an_expr = vpa((subs(Lagrangian_NL, [rho c0 P0 w theta y ], [ 1060 1480 5E4 2*pi*freq0 21 0 ])));
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
syms x z y z t w c0 rho P0 theta beta_mn
clc
freq0 = 1/(an.domain.tpar(2)-an.domain.tpar(1));
gauss_dl                  = 10.0/an.medium.freq0;
gauss_win                 = 5.0/an.medium.freq0;
p_an                =  P0*(sin(w*(t-gauss_dl - sind(theta)*x/c0 - cosd(theta)*z/c0))+ sin(w*(t-gauss_dl + sind(theta)*x/c0 - cosd(theta)*z/c0)));

p_an_expr = vpa((subs(p_an, [rho c0 P0 w y ], [ 1060 1480 1E5 2*pi*freq0 0 ])));
N2 = p_an_expr; N2 = strrep(char(N2),'*','.*');N2  = strrep(char(N2),'^','.^');N2=strrep(char(N2),'/','./')
%%

theta = 25;
c0_an = double(subs(c0,[c0],1480));
T = an.domain.tpar;
X = an.domain.xpar*1E3;
Z = an.domain.zpar*1E3;

len_before = floor(gauss_dl * freq0);
len_window =  floor(gauss_win * freq0);
len_after  = T  - len_window - len_before;
rect_win  = [zeros(len_before,1); rectwin(len_window); zeros(len_after,1)];

[x,t,z]=meshgrid(X,T,Z);
N2 = 100000.0.*sin(471238898.03846895694732666015625.*t - 318404.6608368033492887342298353.*z.*cos(0.017453292519943295769236907684886.*theta) + 318404.6608368033492887342298353.*x.*sin(0.017453292519943295769236907684886.*theta) - 314.15926535897932367815721741733) - 100000.0.*sin(318404.6608368033492887342298353.*z.*cos(0.017453292519943295769236907684886.*theta) - 471238898.03846895694732666015625.*t + 318404.6608368033492887342298353.*x.*sin(0.017453292519943295769236907684886.*theta) + 314.15926535897932367815721741733);
N2 = N2/freq0^2*c0_an^2*.rect_win;
plot_val = 20*log10(squeeze(max(abs(N2))));

figure('WindowState','Maximized');imagesc(Z,X,plot_val);colorbar;
xlabel('Depth [mm]')
ylabel('X [mm]')
title(['Analytical for Source Term due to Lagrangian, @ ', num2str(theta), ' deg'])

xAM_Colormaps(file)
set(gcf,'Color','white')
% set(gca,'FontSize',20)
% xlim([1.94 2.15])
% ylim([-0.2 0.2])

figure('WindowState','Maximized');plot(Z,plot_val(length(X)/2+1,:));
xlabel('Depth [mm]')
ylabel('Pressure')
grid on;
set(gcf,'Color','white')

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
w = 2* pi * freq0 * 3;
theta = 21;

T = an.domain.tpar(30:60);
X = an.domain.xpar*1E3;
Z = an.domain.zpar*1E3;
[x,t,z]=meshgrid(X,T,Z);
blackman_L             = 3;
blackman_delay         = 6; 
blwin                  = 0.35869 - 0.48829 * cos(2*pi*(freq0*t-blackman_delay)/blackman_L) + 0.14128*cos(4*pi*(freq0*t-blackman_delay)/blackman_L) - 0.01168*sin(6*pi*(freq0*t-blackman_delay)/blackman_L);
signature              = (sign(t-blackman_delay/freq0)-sign(t-blackman_delay/freq0-blackman_L/freq0))/2;
an.gauss_dl            = 6.0/an.medium.freq0;
an.gauss_win           = 3.0/an.medium.freq0;
signature              = 1E5*sin(w*(t-gauss_dl)- an.medium.kappa0*cosd(180-theta)*z);
window                 = blwin.*signature ;

N2 = 8*w.^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*t - 2 * w/c0* cosd(theta) * z) + cos(2*w/c0*sind(theta) * x));

% N2 = 8*w.^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*t - 2 * w/c0* cosd(theta) * z) + cos(2*w/c0*sind(theta) * x));
% N2 = 4*w^2*P0^2/rho/c0^4*(sind(theta)^2*(1+cosd(theta)^2)*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-sind(theta)^2*cos(2*w/c0*sind(theta) * x));

N2 = N2/freq0^2*c0^2.*blwin;
N2 = 1E5*(sin(w*(t-gauss_dl) - an.medium.kappa0*cosd(180-theta)*x - an.medium.kappa0*cosd(180-theta)*z) + sin(w*(t-gauss_dl) + an.medium.kappa0*cosd(180-theta)*x - an.medium.kappa0*cosd(180-theta)*z));

plot_val = 20*log10(squeeze(max(abs(N2))));

figure('WindowState','Maximized');imagesc(Z+1,X,plot_val);colorbar;
xlabel('Depth [mm]')
ylabel('X [mm]')
title(['Analytical for Source Term due to Lagrangian, @ ', num2str(theta), ' deg'])
caxis(max(plot_val,[],'all') - [15 0])
xlim([0.8 4])

xAM_Colormaps(file)
set(gcf,'Color','white')
% set(gca,'FontSize',20)
% xlim([1.94 2.15])
% ylim([-0.2 0.2])

figure('WindowState','Maximized');plot(Z,plot_val(length(X)/2+1,:));
xlabel('Depth [mm]')
ylabel('Pressure')
grid on;
set(gcf,'Color','white')

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
theta = 31;

T = an.domain.tpar;
X = an.domain.xpar*1E3;
Z = an.domain.zpar*1E3;
[z,t]=meshgrid(Z,T);
[x,t,z]=meshgrid(X,T,Z);

MediumNL_an = 8*w.^2*P0^2/rho/c0^4*(sind(theta))^2*cos(2*w*t - 2* w/c0 * cosd(theta) * z)-4*w^2*P0^2/rho/c0^4*(sind(theta))^4*(cos(2*w*t - 2 * w/c0* cosd(theta) * z) + cos(2*w/c0*sind(theta) * x));
% MediumNL_an = (beta_mn*((2*P0^2*((11417981541647679048466287755595961091061972992*sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2))/25373292314772621169451365512025 - sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025)^2 + w^2*sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2) + 2*w*cos((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025))^2)/w^2 - (2*P0^2*(sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025) - w*cos((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2))*(sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((32592575621351777380295131014550050576823494298654980010178247189670100796213387298934358016*z*cos((theta*pi)/180))/119103733134816381629145823475159997815219594653880684227626115625 - (260740604970814219042361048116400404614587954389239840081425977517360806369707098391474864128*t)/643803962890899360157544991757621609811997808939895590419600625 + (32592575621351777380295131014550050576823494298654980010178247189670100796213387298934358016*x*sin((theta*pi)/180))/119103733134816381629145823475159997815219594653880684227626115625 + 104296241988325682897342539810539940662930897771211996554195405672790961613448354463744/643803962890899360157544991757621609811997808939895590419600625) + (11417981541647679048466287755595961091061972992*sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025))/25373292314772621169451365512025 - (34253944624943037145398863266787883273185918976*w*cos((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2))/25373292314772621169451365512025 - sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025)^3 - w^3*cos((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2) + 3*w*cos((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025)^2 + 3*w^2*sin((8752240673203471*z*cos((theta*pi)/180))/137438953472 + (8752240673203471*x*sin((theta*pi)/180))/137438953472 - w*(t - 944473296573929/2361183241434822606848))*exp(-((9444732965739290427392*z*cos((theta*pi)/180))/931880319286276675 - (75557863725914323419136*t)/5037190915060955 + (9444732965739290427392*x*sin((theta*pi)/180))/931880319286276675 + 30223145490365728/5037190915060955)^2)*((1427247692705959881058285969449495136382746624*z*cos((theta*pi)/180))/4694059078232934916348502619724625 - (11417981541647679048466287755595961091061972992*t)/25373292314772621169451365512025 + (1427247692705959881058285969449495136382746624*x*sin((theta*pi)/180))/4694059078232934916348502619724625 + 4567192616659071412712425543544147542016/25373292314772621169451365512025)))/w^2))/(c0^4*rho)

MediumNL_an = MediumNL_an/freq0^2*c0^2;

figure('WindowState','Maximized');%plot(Z,squeeze(max(abs(MediumNL_an))))
imagesc(Z,X,squeeze(max(abs(MediumNL_an))))

xlabel('Depth [mm]')
ylabel('Pressure')
grid on;
set(gcf,'Color','white')
% set(gca,'FontSize',20)

%% REMOVE 
theta = 21 
arr.N_el = 32;
arr.H_el = 1;    %mm
arr.W_el = 0.08; %mm
arr.Kerf = 0.02; %mm

arr.W = (arr.W_el + arr.Kerf)* (arr.N_el-1) + arr.W_el;


T = an.domain.tpar;
X = an.domain.xpar*1E3;
Z = an.domain.zpar*1E3;
Z_0 = Z -  Z(1);
X_0 = X -  X(1);
[z,t]=meshgrid(Z,T);
[z,x]=meshgrid(Z_0,X_0);

q1 = x *0;
q2 = q1 ; q3 = q1; q4 = q1; q = q1;

apod_size = X_0(30);
for idxX = 1 : length(X)
    q1(idxX,:) = z(idxX,:)<(X_0(idxX) - apod_size)*tand(90-theta);
    q2(idxX,:) = z(idxX,:)>(X_0(idxX) - max(X_0)/2 + apod_size)*tand(90-theta);
    q3(idxX,:) = z(idxX,:)<(X_0(end) - X_0(idxX) - apod_size)*(tand(90-theta));
    q4(idxX,:) = z(idxX,:)>(X_0(end) - X_0(idxX) - max(X_0)/2 + apod_size)*(tand(90-theta));
    
    q = q1.*q2.*q3.*q4;
end
figure;imagesc(Z,X,q)

