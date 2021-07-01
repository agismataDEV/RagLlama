function aperture = set_ICS_aperture_rectangular(rectwidth,rectheight)

% aperture = set_ICS_aperture_rectangular(rectwidth,rectheight)
% Generates the aperture structure for a rectangular transducer. The
% aperture structure contains the aperture settings used in the 
% set_ICS_source_parametric function. 
% Inputs:
% rectwidth:  width of the transducer in [mm]
% rectheight: height of the transducer in [mm]

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(nargin>2)
    error('Too many arguments')
end

if((rectwidth<=0)||(rectheight<=0))
    error('Invalid value for any of the parameters')
end

aperture = struct('aperturetype','rectangular',...
    'rectwidth',rectwidth,'rectheight',rectheight);