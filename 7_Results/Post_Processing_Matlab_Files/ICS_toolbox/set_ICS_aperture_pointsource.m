function aperture = set_ICS_aperture_pointsource

% aperture = set_ICS_aperture_pointsource 
% Generates the aperture structure for a point source transducer. The
% structure contains the aperture settings used the set_ICS_source_parametric function. 

% JH, version 090116

aperture = struct('aperturetype','pointsource');