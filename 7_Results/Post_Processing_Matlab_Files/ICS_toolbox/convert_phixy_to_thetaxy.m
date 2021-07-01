function thetaxy = convert_phixy_to_thetaxy(phixy)

% thetaxy = convert_INCS_phixy_to_thetaxy(phixy)
% Converts the angles phixy to the angles thetaxy 
% needed in the ICS domainparams structure
% Input:  phixy = [phix phiy]
% phiy:   angle in radians, in the range [-pi/2,pi/2], of the projection of the
%         computational domain on the x-z plane with respect to the positive z-axis, 
%         and positive towards the positive x-axis
% phiy:   angle in radians, in the range [-pi/2,pi/2], of the projection of the
%         computational domain on the y-z plane with respect to the positive z-axis, 
%         and positive towards the positive y-axis
% Output: thetaxy = [thetax thetay]
% thetax: angle in radians, in the range [0 pi], of the projection of the
%         computational domain on the x-z plane, measured with respect to 
%         the positive x-axis
% thetay: angle in radians, in the range [0 pi], of the projection of the
%         computational domain on the y-z plane, measured with respect to
%         the positive y-axis

% JH, version 090116

if(nargin~=1)
    error('Expecting one argument')
end
if numel(phixy)~=2
    error('Expecting a vector of two elements')
end
if(min(phixy)<-pi/2 || max(phixy)>pi/2)
    error('Both elements of phixy should be in the range [-pi/2,pi/2]')
end 

thetaxy = pi/2 - phixy;