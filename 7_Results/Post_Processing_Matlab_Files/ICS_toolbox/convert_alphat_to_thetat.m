function thetat = convert_alphat_to_thetat(alphat)

% thetat = convert_alphat_to_thetat(alphat)
% Converts the factor alphat to the angle thetat
% needed in the ICS domainparams structure
% Input:  
% alphat: factor defining the 'speed' of the comoving 
%         time-frame tau = t - alphat (z/c0)
% Output:
% thetat: angle in radians of the computational domain in the t-z plane
%         with respect to the positive t-axis

% JH, version 090116

if(nargin~=1)
    error('Expecting one argument')
end
if(alphat<0)
    error('Alphat should be zero or positive')
end

if alphat==0
    thetat = pi/2;
else
    thetat = atan(1/alphat);
end