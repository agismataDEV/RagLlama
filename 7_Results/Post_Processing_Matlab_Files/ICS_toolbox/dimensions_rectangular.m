function [phase_delay,traveltime_source,xdim,ydim]=dimensions_rectangular(...
    sourceparams,freq0,c0)

% [phase_delay,traveltime_source,xdim,ydim]=dimensions_rectangular(...
%    sourceparams,freq0,c0)
% Return a number of derived dimensions of the rectangular source aperture
% Inputs:
% sourceparams: structure that has been generated with set_ICS_source_parametric; 
% freq0:        fundamental frequency [Hz]
% c0:           speed of sound in the medium [m/s]
%               it should contain a rectangular aperture.
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
if (~strcmp(sourceparams.aperture.aperturetype,'rectangular'))
    error('Expecting a rectangular source')
end

% shorten notation
ap = sourceparams.aperture;

% Normalize f0 to be counted in MHz, Tpulse to be counted in microseconds, 
% and c0 to be counted in mm/mus
freq0 =  freq0/1e6;
c0 =     c0*1e3/1e6;

%normalize all parameters
rectwidth=ap.rectwidth*freq0/c0;
rectheight=ap.rectheight*freq0/c0;

% obtain the phase delays and travel time over the source
phasedelay = 0;
traveltime_source=sqrt(rectwidth^2+rectheight^2);

% obtain spatial parameters
xdim = rectwidth;
ydim = rectheight;

% set output variables to correct normalization
xdim = xdim * c0/freq0;
ydim = ydim * c0/freq0;

