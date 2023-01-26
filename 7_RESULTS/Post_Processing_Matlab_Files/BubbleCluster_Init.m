% ***********************************************************************************************************
% INITIALIZATION OF THE SIGNIFICANT PARAMETERS
%
% ****************************************************************************
% [medium,domain,slice,file,plin,pnl,Bubble] = BubbleCluster_Init(slicenumber)
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

function [medium,domain,dslice,file,plin,pnl,Bubble] = BubbleCluster_Init(slicenumber)
%SET EVERYTHING TO LATEX INTERPRETER
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% Medium Parameters
medium.freq0               = 1E6;                    %fundamental frequency
medium.c0                  = 1482;                     % Speed of Sound [m/sec]
medium.rho0                = 1060;
medium.dLambdaNN           = medium.c0*1e+3/medium.freq0;
medium.P0                  = 1.01E5  ;

%% Bubble Parameters
Bubble.mindist           = 10e-6;

%% Domain Parameters

domain.beamiterations       = 6;
domain.ibeam                = 0;
domain.PPW_t                = 2;
domain.a_t                  = 0; % Co-moving Time window(1), Non-comoving TIme Window(0)
domain.strides              = [1 1 1 1];
domain.offsets              = [0 0 0 0];
domain.symmetry             = 0;

%% File Parameters
file.rootname            = 'TESTNeumann';
file.dirname             = '../test/1E4P_12mPa_3cycles_1E3MBs_2_4um_9Fnyq';
file.dirname             = '../Paper/Phased_2E5Pa_2micron_1MHz_1E5MBs_90T';
file.dirname             = '../Phased/Phased_1E5MB_1MHz_5E1Pa_2_4um_DHPC';
% file.dirname             = '../ProtonBubbleSync/1E2PS_1E2MB_1_7MHz_1E2Pa_2cycles_regular';
% file.dirname             = '../ProtonBubbleSync/1E0PS_5E3MB_5MHz_1E3Pa_5cycles_regular_mono';
% file.dirname             = '../width_3cycl_round2/PW_2E5Pa_2micron_1MHz_2E4LS';
% file.dirname             = '../width_3cycl_round2/Phased_2E5Pa_2micron_1MHz_1E2MBs';
% file.dirname             = '../../../INCS_XWAVE/7_RESULTS/Microbubbles/David_200deg_x_400kPa_5E5MBs_3_5Fnyq_Z1cm';
% file.dirname             = '../test/Cylinder';
file.contrast_name       = 'ContrastSrc';
file.scatterer           = 'passive_lin';         % 'active','passive_lin', 'passive_nonlin'

file.plot_contrast       = 'yes';            % 'yes' or 'no'
file.plot_attenslices    = 'no';
file.plot_converr        = 'no';
file.plot_colour         = 'viridis';          % 'gray', 'fake_parula' , 'viridis', 'inferno', 'magma', 'plasma'
file.saveplot            = 'no';            % 'yes' or 'no'
file.play_movies         = 'yes';            if (strcmp(file.play_movies,'yes')) ; file.save_movies = 'no'; end
file.load_contrast_from_file ='yes';         % This is to include the Bubble.Contrast inside the BubbleCluster_LocCon source file
file.load_radius_from_file ='no';         % This is to include the Bubble.Contrast inside the BubbleCluster_LocCon source file

%% Slice Parameters
dslice.dims          = [ 'x','y','z'];
dslice.num           = slicenumber;

% Find the files that contain the slice dimension that interests us
LS_slice = split(cellstr(ls([file.dirname,'/',file.rootname int2string_ICS(domain.beamiterations),'_*' int2string_ICS(dslice.num) '_','*.h5'])),[file.rootname int2string_ICS(domain.beamiterations),'_']);
dslice.savedim(dslice.num) = LS_slice{2,end}(1);   % if in cluster change this to LS_slice(2)

% This is the base of the structure of the code
if dslice.savedim(dslice.num) =='y'; domain.dimval = [2 1 3];end
if dslice.savedim(dslice.num) =='x'; domain.dimval = [1 2 3];end
if dslice.savedim(dslice.num) =='z'; domain.dimval = [3 1 2];end

dslice.xlabel = (['\textit{',(dslice.dims(domain.dimval(3))), '} [mm]']);
dslice.ylabel = (['\textit{',(dslice.dims(domain.dimval(2))), '} [mm]']);
dslice.zlabel = (['\textit{',(dslice.dims(domain.dimval(1))), '} [mm]']);

file.focalplanename = [ dslice.savedim(dslice.num) int2string_ICS(dslice.num) ];

%% LINEAR DATA WITHOUT ANY BUBBLES
plin.filename       = [file.dirname '/' file.rootname int2string_ICS(0) '_' file.focalplanename int2string_ICS(domain.ibeam)];

i_start=1;
if isempty(file.dirname); i_start = 2; end

%% LINEAR DATA
output             = load_ICS_slice(plin.filename(i_start:end),domain.strides,domain.offsets);
plin.data          = squeeze(output.data);
domain.tdimpar     = size(plin.data,1);
domain.iProcN      = output.iNumProc;

% %% Source Delay
% file_srcdelay       = [file.dirname '/' file.rootname  'srcdelay' int2string_ICS(0)];
% output_srcdelay     = load_ICS_array(file_srcdelay);
% t_delay             = output_srcdelay.data;
% 
% file_src = [file.dirname '/' file.rootname 'src'];
% output_src      = load_ICS_slice(file_src,domain.strides,domain.offsets);
% psrc       = squeeze(output_src.data);
% domain.tdimpar     = size(psrc,1);
% domain.dimpar      = size(psrc,2);
% domain.zdimpar     = size(psrc,3);
%%
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
domain.Fnyq        = 1/(domain.dtpar*2)/medium.freq0;

%% generate axes for focal plane
domain.TXYstarts=domain.TXYstart+min(domain.TXYoff);

domain.tpar        = (0:domain.tdimpar-1)*domain.dtpar + domain.tstart*domain.dtpar;
domain.par{1}      = (domain.start(1)+(0:domain.dimlen(1)-1))*domain.dxpar;
domain.par{2}      = (domain.start(2)+(0:domain.dimlen(2)-1))*domain.dypar;
domain.par{3}      = (domain.start(3)+(0:domain.dimlen(3)-1))*domain.dzpar;

dslice.pos(dslice.num) = domain.par{domain.dimval(1)}(dslice.index(dslice.num)) ;

%% NONLINEAR DATA
pnl.filename       = [file.dirname '/' file.rootname int2string_ICS(domain.beamiterations) '_' file.focalplanename int2string_ICS(domain.ibeam)];

output_nonl      = load_ICS_slice(pnl.filename);
pnl.data       = squeeze(output_nonl.data);

%% CONTRAST DATA WITHOUT ANY BUBBLES
if (strcmp(file.plot_contrast,'yes') && domain.beamiterations>0)
    pcon.filename       = [file.dirname '/' file.contrast_name int2string_ICS(domain.beamiterations) '_' file.focalplanename int2string_ICS(domain.ibeam)];
    i_start=1;
    if isempty(file.dirname) ;i_start = 2; end
    
    output      = load_ICS_slice(pcon.filename(i_start:end),domain.strides,domain.offsets);
    pnl.contrastdata        = squeeze(output.data);
    
end
if (domain.symmetry~=0) 
    domain.par{domain.symmetry} = [-flip(domain.par{domain.dimval(domain.symmetry)}(2:end)) domain.par{domain.dimval(domain.symmetry)}]; 
    plin.data  = [flip(plin.data(:,(2:end),:) ,domain.symmetry)  plin.data];
    pnl.data  = [flip(pnl.data(:,(2:end),:) ,domain.symmetry)  pnl.data];
end

%% CREATE SAVE FOLDER
file.savedir             = ['../Images/meeting_',num2str(yyyymmdd(datetime('today'))),'/',dslice.savedim(dslice.num),'_',num2str(dslice.pos(dslice.num))];
if strcmp(file.saveplot,'yes')
    mkdir(file.savedir)
end
end