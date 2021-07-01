function [phase_delay,traveltime_source,xdim,ydim]=dimensions_phasedarray(...
    sourceparams,freq0,c0)

% [phase_delay,traveltime_source,xdim,ydim]=dimensions_phasedarray(sourceparams,freq0,c0)
% Return a number of derived dimensions of the phased array source aperture
% Inputs:
% sourceparams: structure that has been generated with set_ICS_source_parametric; 
% freq0:        fundamental frequency [Hz]
% c0:           speed of sound in the medium [m/s]
%               it should contain a phased array aperture.
% Outputs:
% phase_delay:  maximum phase delay on any of the elements in [1/freq0]
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
if (~strcmp(sourceparams.aperture.aperturetype,'phasedarray'))
    error('Expecting a phased array source')
end

% shorten notation
ap = sourceparams.aperture;

% Normalize f0 to be counted in MHz, Tpulse to be counted in microseconds, 
% and c0 to be counted in mm/mus
freq0 =  freq0/1e6;
c0 =     c0*1e3/1e6;

%normalize all parameters to wavelengths
elwidth=ap.elwidth*freq0/c0;
elheight=ap.elheight*freq0/c0;
kerf=ap.kerf*freq0/c0;
focusx=ap.focusx*freq0/c0;
focusz=ap.focusz*freq0/c0;
elevationfocusz=ap.elevationfocusz*freq0/c0;

% obtain the phase delays and travel time over the source
half_arraywidth=.5*(ap.numels*(elwidth+kerf)-kerf);
focusx=abs(focusx);
%phase delay of the leftmost element
max_phase_delay_x=sqrt((focusx+half_arraywidth)^2+focusz^2)-focusz;
if (focusx>half_arraywidth)
    %subtract the phase delay of the rightmost element if
    %focusx>half_arraywidth
    max_phase_delay_x=max_phase_delay_x - ...
        (sqrt((focusx-half_arraywidth)^2+focusz^2)-focusz);
end
if(elevationfocusz>1e-5)
    max_phase_delay_y=sqrt((elheight/2)^2+elevationfocusz^2)-elevationfocusz;
else
    max_phase_delay_y=0;
end

% set output variables
phase_delay = max_phase_delay_x+max_phase_delay_y;
traveltime_source=sqrt((2*half_arraywidth)^2+elheight^2);
xdim = 2*half_arraywidth;
ydim = elheight;

% set output variables to correct normalization
xdim = xdim * c0/freq0;
ydim = ydim * c0/freq0;

