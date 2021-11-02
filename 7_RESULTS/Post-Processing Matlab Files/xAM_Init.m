% ***********************************************************************************************************
% INITIALIZATION OF THE SIGNIFICANT PARAMETERS
%
% ****************************************************************************
% [medium,domain,slice,file,plin,pnl,Bubble] = xAM_Init(slicenumber)
% ****************************************************************************
% Initializes the parameters for the medium, domain, slice, file ,linear pressure field
% non-linear pressure field and the scatterers (Bubbles) structures.
% Loads the files for the linear pressure field and the pressure field created
% from the scatterers and assign the important parameters to the related structure
%
% The visualization code is based on the idea of the dimension of the dslice.
% So it can be used independent of the dimension of the dslice.
% if the slice is from y , then code dimensions become ['y' 'x' 'z'] = [2 1 3]
% if the slice is from x , then code dimensions become ['x' 'y' 'z'] = [1 2 3]
% if the slice is from z , then code dimensions become ['z' 'x' 'y'] = [3 1 2]
%
% Inputs:
% slicenumber   integer     Number of the slice of the domain that contains all the information
%
% Outputs:
% medium        struct      Contains all the medium properties, e.g Density, Speed of Sound etc
% domain        struct      Contains all the computational domain properties, e.g space and time step
% slice         struct      Contains all the information about the slice ,e.g Dimension, Position, Index etc
% file          struct      Contains all the loaded file properties ,e.g directory and root name etc
% plin          struct      Contains all the pressure values of the plin slice
% pnl           struct      Contains all the pressure values of the pnl slice
% Bubble        struct      Contains all properties of the Scatterers' Cluster
%
% ******************
% AM, version 200723
% ******************
%
% ***********************************************************************************************************

function [medium,domain,dslice,file,plin,pnl,Bubble] = xAM_Init(slicenumber)
%% Medium Parameters
medium.freq0               = 15E6;                    %fundamental frequency
medium.c0                  = 1480;                     % Speed of Sound [m/sec]
medium.rho0                = 1000;
medium.dLambdaNN           = medium.c0*1e+3/medium.freq0;
medium.P0                  = 1.01E5  ;

%% Bubble Parameters
Bubble.N                 = 0 ;
Bubble.mindist           = 10e-6;

%% Domain Parameters

domain.beamiterations       = 1;
domain.ibeam                = 0;
domain.Fnyq                 = 3.5;
domain.PPW_t                = 2;
domain.strides              = [1 1 1 1];
domain.offsets              = [0 0 0 0];
domain.symmetry             = 0;

%% File Parameters
file.lin_name            = 'TESTNeumann';
file.dirname             = '../David_210deg_x_Lagrangian_Only';
file.dirname_left        = '../David_150deg_left_high_amp';
file.dirname_right       = '../David_150deg_left_high_amp';
file.rootname            = 'TESTNeumann';
file.scatterer           = 'passive_nonlin';         % 'active','passive_lin', 'passive_nonlin'

file.plot_contrast       = 'no';            % 'yes' or 'no'
file.plot_attenslices    = 'no';
file.plot_converr        = 'yes';
file.plot_colour         = 'viridis';        % 'gray', 'fake_parula' , 'viridis', 'inferno', 'magma', 'plasma'
file.saveplot            = 'no';            % 'yes' or 'no'
file.play_movies         = 'no';            if (strcmp(file.play_movies,'yes')) ; file.save_movies = 'no'; end
file.load_contrast_from_file ='no';         % This is to include the Bubble.Contrast inside the BubbleCluster_LocCon source file

domain.angle = str2num(file.dirname(strfind(file.dirname,'deg')+[-3:-1]))/10;
%% Slice Parameters
dslice.dims          = [ 'x','y','z'];
dslice.num           = slicenumber;

% Find the files that contain the slice dimension that interests us
LS_slice = split(cellstr(ls([file.dirname,'/',file.lin_name int2string_ICS(domain.beamiterations),'_*' int2string_ICS(dslice.num) '_','*.h5'])),[file.lin_name int2string_ICS(domain.beamiterations),'_']);
dslice.savedim(dslice.num) = LS_slice{2,2}(1);   % if in cluster change this to LS_slice(2)

% This is the base of the structure of the code
if dslice.savedim(dslice.num) =='y'; domain.dimval = [2 1 3];end
if dslice.savedim(dslice.num) =='x'; domain.dimval = [1 2 3];end
if dslice.savedim(dslice.num) =='z'; domain.dimval = [3 1 2];end

dslice.xlabel = ([upper(dslice.dims(domain.dimval(3))) ' [mm]']);
dslice.ylabel = ([upper(dslice.dims(domain.dimval(2))) ' [mm]']);
dslice.zlabel = ([upper(dslice.savedim(dslice.num)) ' [mm]']);

%%
% Find the all the log files and split them based on file.rootname to get all possibles processor files .
% Split with the extension and get the max value added +1 because the first processor is 0 , to get the ProcNum
% LS_log = split(cellstr(ls([file.dirname,'/*.log'])),[file.lin_name,'_']);
% LS_logsplit = split(LS_log(2:end,2),'.log');  % if in cluster change this to LS_log(2:end)
% domain.iProcN              = max(str2double(LS_logsplit(:,1)))+1;

file.focalplanename = [ dslice.savedim(dslice.num) int2string_ICS(dslice.num) ];

%% X_wave
plin.filename       = [file.dirname '/' file.lin_name int2string_ICS(0) '_' file.focalplanename int2string_ICS(domain.ibeam)];

i_start=1;
if isempty(file.dirname); i_start = 2; end

%% LINEAR DATA
output             = load_ICS_slice(plin.filename(i_start:end),domain.strides,domain.offsets);
plin.data          = squeeze(output.data);
domain.tdimpar     = size(plin.data,1);

%% Source Delay
% file_srcdelay       = [file.dirname '/' file.rootname  'srcdelay' int2string_ICS(0)];
% output_srcdelay     = load_ICS_array(file_srcdelay);
% t_delay             = output_srcdelay.data;

%%
% domain.iProcN = output.iNumProc;
domain.dimlen(domain.dimval(1)) = output.lengthslice;
domain.dimlen(domain.dimval(2)) = size(plin.data,2);
domain.dimlen(domain.dimval(3)) = size(plin.data,3);

domain.TXYoff      = output.TXYoff; if dslice.savedim(dslice.num) =='z' ; domain.TXYoff = repmat(domain.TXYoff,domain.dimlen(3),1); end
domain.TXYstart    = output.start(1:3);
dslice.index       = output.sliceindex+1;
domain.dtpar       = output.stepsize(1);
domain.dxpar       = output.stepsize(2)*1e+3; % step size dx in mm
domain.dypar       = output.stepsize(3)*1e+3; % step size dy in mm
domain.dzpar       = output.stepsize(4)*1e+3; % step size dz in mm
domain.start(1)    = domain.TXYstart(2)+domain.TXYoff(1,2);
domain.start(2)    = domain.TXYstart(3)+domain.TXYoff(1,3);
domain.start(3)    = output.start(4);
domain.tstart      = domain.TXYstart(1)+domain.TXYoff(1,1);
domain.xsource     = (output.start(2)+(0:size(plin.data,2)-1))*domain.dxpar;
domain.ysource     = (output.start(3)+(0:size(plin.data,3)-1))*domain.dypar;

%% generate axes for focal plane
domain.TXYstarts=domain.TXYstart+min(domain.TXYoff);

domain.tpar        = (0:domain.tdimpar-1)*domain.dtpar;
domain.par{1}      = (domain.start(1)+(0:domain.dimlen(1)-1))*domain.dxpar;
domain.par{2}      = (domain.start(2)+(0:domain.dimlen(2)-1))*domain.dypar;
domain.par{3}      = (domain.start(3)+(0:domain.dimlen(3)-1))*domain.dzpar;

dslice.pos(dslice.num) = domain.par{domain.dimval(1)}(dslice.index(dslice.num)) ;

%% NONLINEAR DATA
pnl.filename       = [file.dirname '/' file.rootname int2string_ICS(domain.beamiterations) '_' file.focalplanename int2string_ICS(domain.ibeam)];

output_nonl      = load_ICS_slice(pnl.filename);
pnl.xwave_data       = squeeze(output_nonl.data);

%% CONTRAST DATA WITHOUT ANY BUBBLES
if strcmp(file.plot_contrast,'yes')
    pcon.filename       = [file.dirname_right,'/', file.rootname int2string_ICS(domain.beamiterations) '_' file.focalplanename int2string_ICS(domain.ibeam)];    
    output      = load_ICS_slice(pcon.filename(i_start:end),domain.strides,domain.offsets);
    pnl.rightwave_data        = squeeze(output.data);
   
    pcon.filename       = [file.dirname_left,'/', file.rootname int2string_ICS(domain.beamiterations) '_' file.focalplanename int2string_ICS(domain.ibeam)];
    output      = load_ICS_slice(pcon.filename(i_start:end),domain.strides,domain.offsets);
    pnl.leftwave_data        = squeeze(output.data);
    
end
if (domain.symmetry~=0) 
    domain.par{domain.symmetry} = [-flip(domain.par{domain.dimval(domain.symmetry)}(2:end)) domain.par{domain.dimval(domain.symmetry)}]; 
    plin.data  = [flip(plin.data(:,(2:end),:) ,domain.symmetry)  plin.data];
    pnl.xwave_data  = [flip(pnl.xwave_data(:,(2:end),:) ,domain.symmetry)  pnl.xwave_data];
end

%% CREATE SAVE FOLDER
file.savedir             = ['../Images/meeting_',num2str(yyyymmdd(datetime('today'))),'/',dslice.savedim(dslice.num),'_',num2str(dslice.pos(dslice.num))];
if strcmp(file.saveplot,'yes')
    mkdir(file.savedir)
end
end