function aperture = set_ICS_aperture_triangular(triangedge,triangxcenter)

% aperture = set_ICS_aperture_triangular(triangedge,triangxcenter)
% Generates the aperture structure for a triangular transducer. 
% The top corner of the triangle is positioned on the positive x-axis.
% The aperture structure contains the aperture settings used in
% the set_ICS_source_parametric function. 
% Inputs:
% triangedge:    length of each of the three edges in [mm]
% triangxcenter: position of the center on the x axis in [mm]

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(nargin>2)
    error('Too many arguments')
end

if((triangedge<=0))
    error('Invalid value for triangedge')
end

aperture = struct('aperturetype','triangular',...
    'triangedge',triangedge,'triangxcenter',triangxcenter);