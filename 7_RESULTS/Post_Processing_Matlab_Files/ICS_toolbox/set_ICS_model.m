function modelparams = set_ICS_model(freq0,Fnyq,beamiterations,...
    contrastsourcetype,slicesavetype,slicesavedim,slicesavepos,...
    debuglevel,contrastsourceflags,FDorders)

% modelparams = set_ICS_model(freq0,Fnyq,beamiterations,...
%     contrastsourcetype,slicesavetype,slicesavedim,slicesavepos,...
%     debuglevel,contrastsourceflags,FDorders)
% Generates the modelparams structure that contains all general model settings
% used in the ICS Parnac program. 
% Inputs:
% freq0:               Fundamental frequency [Hz]
% Fnyq:                Discretization frequency factor [relative to freq0]
% beamiterations:      vector with the number of iterations for each beam;
%                      all numbers should be 1 or larger, and not increasing
% contrastsourcetype:  (optional, default = 'nonlin') Type of the contrast source, 
%                      choose between 'nonlin', 'complexcontrast', 'sphere', 
%                      'luneberg' or 'blob'
% slicesavetype:       (optional, default = 'firstandlast') of which iterations the
%                      slices are saved, choose between 'firstandlast', 'last' or 'all'
% slicesavedim:        (optional, default = 'y') vector with the dimension for each 
%                      slice to be saved, choose between 'x', 'y' or 'z'
% slicesavepos:        (optional, default = 0) vector with the positions in mm 
%                      in the dimension slicesavedim(:) for each slice to be saved
% debuglevel:          (optional, default = 4) Level of feedback through the log files 
% contrastsourceflags: (optional, default = [true true true]) Three flags 
%                      [useantialiasing usesupporttapering usefreqtapering]
%                      for the evaluation of the contrast source operator
% FDorders:            (optional, default = [30 30]) Overall order of the 
%                      Finite Difference schemes in the temporal and spatial
%                      dimensions -- should be even numbers

% JH, version 090116

if(nargin<3)
    error('Not enough arguments')
end
if(nargin<4 || isempty(contrastsourcetype))
    contrastsourcetype = 'nonlin';
end
if(nargin<5 || isempty(slicesavetype))
    slicesavetype = 'firstandlast';
end
if(nargin<6 || isempty(slicesavedim))
    slicesavedim = 'y';
end
if(nargin<7 || isempty(slicesavepos))
    slicesavepos = 0;
end
if(nargin<8 || isempty(debuglevel))
    debuglevel = 4;
end
if(nargin<9 || isempty(contrastsourceflags))
    contrastsourceflags = [true true true];
end
if(numel(contrastsourceflags)~=3)
    error('contrastsourceflags should be a vector of three elements')
end
if(nargin<10 || isempty(FDorders))
    FDorders = [30 30];
end
if(numel(FDorders)~=2)
    error('FDorders should be a vector of two elements')
end
if(nargin>10)
    error('Too many arguments')
end

if((freq0<=0)||(Fnyq<1))
    error('freq0 of Fnyq have invalid values')
end
if(isempty(strmatch(contrastsourcetype,strvcat(...
        'nonlin','complexcontrast','sphere','luneberg','blob'), 'exact')))
    error('contrastsourcetype has an invalid value')
end
if((numel(beamiterations)<1)||(min(beamiterations)<1))
    error('beamiterations should be a vector of numbers >=1')
end
if(isempty(strmatch(slicesavetype,strvcat(...
        'firstandlast','last','all'), 'exact')))
    error('slicesavetype has an invalid value')
end

if(numel(slicesavedim)~=numel(slicesavepos))
    error('slicesavedim should have the same number of elements as slicesavepos')
end
for i=1:numel(slicesavedim)
    if(isempty(strmatch(slicesavedim(i),strvcat('x','y','z'),'exact')))
        error('Elements of slicesavedim should be ''x'',''y'' or ''z''')
    end
end

debuglevel=round(debuglevel);

if(~islogical(contrastsourceflags))
    error('contrastsourceflags should contain logical values, i.e. false or true')
end

modelparams = struct('freq0',freq0,'Fnyq',Fnyq,...
    'beamiterations',beamiterations,'contrastsourcetype',contrastsourcetype,...
    'slicesavetype',slicesavetype,'slicesavedim',slicesavedim,...
    'slicesavepos',slicesavepos,'debuglevel',debuglevel,...
    'useantialiasing',contrastsourceflags(3),...
    'usesupporttapering',contrastsourceflags(2),...
    'usefreqtapering',contrastsourceflags(3),...
    'FDTorder',FDorders(1),'FDXorder',FDorders(2));