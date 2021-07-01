function mediumparams = set_ICS_medium(rho0,c0,beta,a_att,b_att,fc0)

% mediumparams = set_ICS_medium(p0,rho0,c0,beta,a_att,b_att,fc0)
% Generates the mediumparams structure that contains all information on the
% medium used in the ICS Parnac program.
% Inputs:
% rho0:  Ambient mass density [kg/m^3]
% c0:    Small signal sound velocity [m/s]
% beta:  (optional, default = 0) Nonlinearity parameter
% a_att: (optional, default = 0) Frequency-dependent attenuation parameter [Np/cm MHz^b]
% b_att: (optional, default = 2) Power in frequency power attenuation law 
% fc0:   (optional, default = -1) Frequency at which c0 is given [Hz] -- if 
%                                 fc0 = -1, then it is set to f0 in
%                                 generate_ICS_configurationfile()

%        
% JH, version 090116

if(nargin<2)
    error('Not enough arguments')
elseif(nargin<3 || isempty(beta))
    beta = 0;
end
if(nargin<4 || isempty(a_att))
    a_att = 0;
end
if(nargin<5 || isempty(b_att))
    b_att = 2;
end
if(nargin<6 || isempty(fc0))
    fc0 = -1;
end
if(nargin>6)
    error('Too many arguments')
end

if((rho0<=0)||(c0<=0)||(beta<0)...
        ||(a_att<0)||(b_att<=1)||(b_att>2)||((fc0~=-1)&&(fc0<0)))
    error('One of the parameters has a negative or otherwise invalid value')
end

mediumparams = struct('rho0',rho0, 'c0', c0, 'beta', beta, ...
    'a_att', a_att, 'b_att', b_att, 'fc0', fc0);