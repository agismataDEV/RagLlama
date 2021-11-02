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
    phi_lagr = 1/rho * int(p_lagr,t);
    phi_sq = simplify(phi_lagr^2);
    phi_sq_x = simplify(diff(diff(phi_sq,x),x));
    phi_sq_y = simplify(diff(diff(phi_sq,y),y));
    phi_sq_z = simplify(diff(diff(phi_sq,z),z));
    phi_sq_t = simplify(diff(diff(phi_sq,t),t));
    Lagrangian_simpl = subs(rho*simplify(phi_sq_x + phi_sq_y + phi_sq_z - 1/c0^2 * phi_sq_t)/4 , theta, 0  );
    L_String = strjoin({char("The pressure due to Lagrangian equals to : "),char(Lagrangian_simpl)});
    disp(L_String)
    
    phi_sq_x = simplify(diff(phi_lagr,x)^2);
    phi_sq_y = simplify(diff(phi_lagr,y)^2);
    phi_sq_z = simplify(diff(phi_lagr,z)^2);
    phi_sq_t = simplify(diff(phi_lagr/c0,t)^2);
    Lagrangian = rho*simplify(phi_sq_x + phi_sq_y + phi_sq_z - phi_sq_t)/2;
    Lagrangian_sub = subs(Lagrangian , theta, 21  );
    Lagrangian_sub = subs(Lagrangian_sub , [rho c0 P0 w theta x z ], [1060 1480 4E4 2*pi*freq0 21 0 domain.parz(143)*1E-3]  );
    L_String = strjoin({char("The Lagrangian equals to : "),char(Lagrangian_sub)});
    disp(L_String)
    
    L1_X = diff(diff(Lagrangian,x),x);
    L1_Y = diff(diff(Lagrangian,y),y);
    L1_Z = diff(diff(Lagrangian,z),z);
    L2_T = diff(diff(Lagrangian,t),t)/c0^2;
    Lagrangian_NL = L1_X+L1_Y+L1_Z+L2_T;
    L_String = strjoin({char("The nonlinear term due to Lagrangian equals to : "),char(Lagrangian_NL)});
%     Lagrangian_NL_sub = vpa(simplify(subs(Lagrangian_NL , [P0 theta c0 rho x y z w ],[1E6 15 1480 1060 0 0 domain.par{3}(300)*1E-3 2*pi*8*15E6] )));
%     Lagrangian_NL_sub_zt = subs(Lagrangian_NL_sub, t ,domain.tpar);
%     L_String = strjoin({char("The nonlinear term due to Lagrangian equals to : "),char(Lagrangian_NL_sub)});
figure;plot(domain.part+domain.parz(142)*1E-3/1480+max(td_x)*1E-6,real(double(subs(Lagrangian_sub * phi.blackman_win * phi.sign_win,t,domain.part))))
    disp(L_String)
    
    Medium_NL = diff(diff(p_lagr^2,t),t)/c0^4/rho*beta_mn;
    L_String = strjoin({char("The nonlinear term due to Medium nonlinearities equals to : "),char(Medium_NL)});
%     Medium_NL_sub = vpa(simplify(subs(Medium_NL , [P0 theta c0 rho x y z w beta_mn ],[1E6 15 1480 1060 0 0 domain.par{3}(300)*1E-3 2*pi*8*15E6 3] )));
%     Medium_NL_sub_zt = subs(Medium_NL_sub, t ,domain.tpar);
%     L_String = strjoin({char("The nonlinear term due to Medium nonlinearities equals to : "),char(Medium_NL_sub)});
    disp(L_String)
end
