function aperture = set_ICS_aperture_cylindrical(radius,radfocus)

% aperture = set_ICS_aperture_cylindrical(radius,radfocus)
% Generates the aperture structure for a cylindrical transducer. The
% aperture structure contains the aperture settings used in the 
% set_ICS_source_parametric function. 
% Inputs:
% radius:   radius of the transducer in [mm]
% radfocus: (optional, default = 0) Focal distance in mm for a focused
%           or defocused cylindrical transducer. Zero means unfocused.

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(nargin<2 || isempty(radfocus))
    radfocus = 0;
end
if(nargin>2)
    error('Too many arguments')
end

if(radius<=0)
    error('Invalid value for the radius')
end

aperture = struct('aperturetype','cylindrical',...
    'radius',radius,'radfocus',radfocus);