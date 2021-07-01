% *********************************************************************************
% MARMOTTANT MODEL SOLVER
%
% **********************************************************
% [dR] = RP(t,R,t_p_driv,p_driv,par)
% **********************************************************
% Solves the Rayleigh - Plesset Equation with the addition of the
% Marmottant Model
%
% Inputs:
% t          real      Time stamp for the solution of RP model
% R          real      Array with initial values for the solution
% p_driv     real      Array with the values of the driving pressure
% t_p_driv   real      Array with the values of the temporal dimension of p_driv
%
% Outputs:
% dR         real      Solution of RP model .
%                      Contains the values of the all the time derivatives
%
% ******************
% AM, version 200723
% ******************
%
% *********************************************************************************

function [dR] = RP(t,R,t_p_driv,p_driv,par)
% R(1) = R
% R(2) = R_dot

p_driv = interp1(t_p_driv,p_driv,t);  %Find pressure of the driving pulse at time t when p_driv(t_p_driv) is known

%%======================= Marmotant model - Shell buckling and rupture ====

if R(1)<par.R_b                                %buckled state
    par.sigma_R = 0;
elseif R(1)>par.R_r                            %rupture state
    par.sigma_R = par.sigma_w ;
else                                           %Elastic state
    par.sigma_R = par.chi*(R(1)^2/par.R_b^2-1) ;
end

P_elas =  2*par.sigma_R/R(1) ;
P_vis  =  4*par.S_vis*R(2)/R(1)^2 ;

%%===============================Solving for R_dd ========================

dR = zeros(2,1);           % Create column vector for the output
dR(1)= R(2) ;              % dR(1) = R_dot
dR(2)=((par.P_g0*(par.R0/R(1))^(3*par.gamma)*(1-3*par.gamma*R(2)/par.c)-par.P0-p_driv-4*par.mu*R(2)/R(1)-P_elas-P_vis)/par.rho-3/2*R(2)^2)/R(1) ;
