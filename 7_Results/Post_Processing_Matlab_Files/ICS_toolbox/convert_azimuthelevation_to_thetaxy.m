function thetaxy = convert_azimuthelevation_to_thetaxy(azimuthelevation)

% thetaxy = convert_azimuthelevation_to_thetaxy(azimuthelevation)
% Converts azimuth and elevation to the angles thetaxy 
% needed in the ICS domainparams structure
% Input: azimuthelevation = [azimuth, elevation]
% azimuth:   rotation angle in radians in the range [-pi/2,pi/2] of the 
%            main axis of the computational domain around the y-axis. 
%            An azimuth of 0 radians locates the main axis
%            in the y-z plane, and a positive azimuth turns it towards
%            the positive x-axis.
% elevation: tilt angle in radians in the range [-pi/2,pi/2] of the 
%            computational domain off the x-z plane, positive in the
%            direction of the positive y-axis.
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
if numel(azimuthelevation)~=2
    error('Expecting a vector of two elements')
end

azimuth = azimuthelevation(1);
elevation = azimuthelevation(2);

if(azimuth<-pi/2 || azimuth>pi/2)
    error('Azimuth should be in the range [-pi/2 pi/2]')
end
if(elevation<-pi/2 || elevation>pi/2)
    error('Elevation should be in the range [-pi/2 pi/2]')
end

thetaxy(1) = pi/2 - azimuth;
if(abs(azimuth)==pi/2)
    thetaxy(2)=0;
else
    thetaxy(2) = pi/2 - atan(tan(elevation)/cos(azimuth));
end
