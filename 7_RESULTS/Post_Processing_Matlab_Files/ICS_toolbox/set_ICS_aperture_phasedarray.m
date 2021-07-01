function aperture = set_ICS_aperture_phasedarray(numels,elwidth,elheight,kerf,...
    focusx,focusz,elevationfocusz,apodization,apofilename)

% aperture = set_ICS_aperture_phasedarray(numels,elwidth,elheight,kerf,...
%     focusx,focusz,elevationfocusz,apodization) 
% Generates the aperture structure for a phased array transducer, and if applicable
% the apodization file. The aperture structure contains the aperture 
% settings used in the set_ICS_source_parametric function. 
% Inputs:
% numels:         number of elements
% elwidth:        width of each element [mm]
% elheight:       height of each element [mm]
% kerf:           distance between the elements [mm]
% focusx:         x-coordinate of the focus [mm]
% focusz:         z-coordinate of the focus [mm]
% elevationfocus: z-coordinate of the focus in the elevation focusing,
%                 i.e. focusing in the y-direction
% apodization:    (optional, default = ones for all elements) vector 
%                 of numel elements containing the apodization of the elements
% apofilename:    (optional, default = 'apodization') filename of the 
%                 apodization file, without the extension.

% JH, version 090116

if(nargin<7)
    error('Not enough arguments')
end
if(nargin<8 || isempty(apodization))
    apodization = ones(1,numels);
end
if(nargin<9 || isempty(apofilename))
    apofilename = 'apodization';
end
if(nargin>9)
    error('Too many arguments')
end

if((numels<1)||(elwidth<=0)||(elheight<=0)||(kerf<0)||(focusz<=0)...
    ||(elevationfocusz<=0))
    error('Invalid value for any of the parameters')
end
if(numel(apodization)~=numels)
    error('Number of elements of the apodization vector should equal numels')
end

if(any(apodization~=1))

    apofilename = strcat(apofilename,'.dat');
    sizeapo=size(apodization);
    save(apofilename,'-ASCII','-DOUBLE','sizeapo','apodization');

else

    apofilename = 'none';

end

aperture = struct('aperturetype','phasedarray',...
    'numels',numels,'elwidth',elwidth,'elheight',elheight,'kerf',kerf,...
    'focusx',focusx,'focusz',focusz,'elevationfocusz',elevationfocusz,...
    'apodization',apodization,'apofilename',apofilename);