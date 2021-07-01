function [phase_delay,traveltime_source,xdim,ydim]=dimensions_cylindrical(...
    sourceparams,freq0,c0)

% [phase_delay,traveltime_source,xdim,ydim]=dimensions_cylindrical(...
%    sourceparams,freq0,c0)
% Return a number of derived dimensions of the cylindrical source aperture
% Inputs:
% sourceparams: structure that has been generated with set_ICS_source_parametric; 
% freq0:        fundamental frequency [Hz]
% c0:           speed of sound in the medium [m/s]
%               it should contain a cylindrical aperture.
% Outputs:
% phase_delay:  maximum phase delay on the aperture in [1/freq0]
% traveltime_source: travel time needed to cross the source at its maximum
%               extension in [1/freq0]
% xdim:         dimension in x in [mm]
% ydim:         dimension in y in [mm]

if(nargin~=3)
    error('Expecting three arguments')
end
if (~strcmp(sourceparams.definitiontype,'parametric'))
    error('Expecting a parametric source definition')
end
if (~strcmp(sourceparams.aperture.aperturetype,'cylindrical'))
    error('Expecting a cylindrical source')
end

% shorten notation
ap = sourceparams.aperture;

% Normalize f0 to be counted in MHz, Tpulse to be counted in microseconds, 
% and c0 to be counted in mm/mus
freq0 =  freq0/1e6;
c0 =     c0*1e3/1e6;

%normalize all parameters
radius=ap.radius*freq0/c0;
radfocus=ap.radfocus*freq0/c0;

% obtain the phase delays and travel time over the source
if (abs(radfocus)>1e-5)
    phase_delay=sqrt(radius^2+radfocus^2)-abs(radfocus);
else
    phase_delay=0;
end

traveltime_source=2*radius;

% obtain spatial parameters
xdim = 2*radius;
ydim = 2*radius;

% set output variables to correct normalization
xdim = xdim * c0/freq0;
ydim = ydim * c0/freq0;

