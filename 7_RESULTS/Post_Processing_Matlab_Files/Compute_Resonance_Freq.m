c = 1480;
chi  = 0.7;
mu = 2E-3;
rho = 1060;
gamma = 1.07;
sigma_R0 = 0.00;
kappa_s = 3E-8;
R  = 2E-6;
f0 = 1/(2*pi*R*sqrt(rho))*sqrt(3*gamma*1.01E5 + (3*gamma-1)*2*sigma_R0/R  + 4*chi/R)

omega = 3*1.01E5*gamma/R^2/rho+4*chi/R^3/rho;
delta = R/c + 4*mu/(rho*R^2*omega) + 4*kappa_s/(rho*R^3*omega);

fr = f0*sqrt(1-delta^2/2)