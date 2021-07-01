function [memsize,numpoints] = check_ICS_memory(mediumparams,domainparams,...
    modelparams,sourceparams,numprocs,showoutput,dTdelay,iDimSource)

% [memsize,numpoints] = check_ICS_memory(mediumparams,domainparams,...
%    modelparams,sourceparams,numprocs,showoutput,dTdelay,iDimSource)
% Calculates the memorysizes per processor for several critical instances in the
% ICS Parnac program. The function emulates the Fortran code that determines the 
% arrays and array sizes in the program.
% Inputs:
% mediumparams: structure that has been generated with set_ICS_medium.
% domainparams: structure that has been generated with set_ICS_domain.
% modelparams:  structure that has been generated with set_ICS_model.
% sourceparams: structure that has been generated with set_ICS_source_parametric
%               or with set_ICS_source_data.
% numprocs:     number of processors
% showoutput:   (optional, default = 0) show detailed summary of number of
%               points (1) or not (0)
% dTdelay:      (optional, default = 0) dTdelay of the source data input file
%               if any is used
% iDimSource:   (optional, default = [0 0 0]) vector of three numbers giving
%               the T/X/Y dimension sizes of the source data input file.
%               If a source signature input file is used, only the T
%               dimension is used, the rest is derived from the parametric
%               aperture.
% Outputs:
% memsize:      vector with the critical memorysizes per processor in MB
% numpoints:    11x4 element array containing the number of points in the
%               T/X/Y/Z dimensions. The rows correspond with the
%               following arrays employed internally in the ICS program:
%               1 - The field beam
%               2 - The source slice
%               3 - The field  beam including d/dz margins for source2beam
%               4 - The plane slice
%               5 - The field  beam including d/dz margins for plane2beam
%               6 - The source beam slice for the beam-to-beam calculations
%               7 - The field  beam slice for the beam-to-beam calculations
%               8 - The source beam slice for the source-slice-to-beam calculations
%               9 - The field  beam slice for the source-slice-to-beam calculations
%               10 - The source beam slice for the plane-to-beam calculations
%               11 - The field  beam slice for the plane-to-beam calculations

if(nargin<3)
    error('Not enough arguments')
end
if(nargin<6 || isempty(showoutput))
    showoutput = 0;
end
if(nargin<7 || isempty(dTdelay))
    dTdelay = 0;
end
if(nargin<8 || isempty(iDimSource))
    iDimSource = [0 0 0];
end
if(nargin>8)
    error('Too many arguments')
end

np = get_ICS_arraysizes(mediumparams,domainparams,...
    modelparams,sourceparams,numprocs,dTdelay,iDimSource);

if(showoutput ~=0)
    disp('Number of points in T/X/Y/Z in the field beam:')
    disp([' ' num2str(np(1,:))]);
    disp('Number of points in T/X/Y/Z in the source slice:')
    disp([' ' num2str(np(2,:))]);
    disp('Number of points in T/X/Y/Z in the field beam with d/dz margins:')
    disp([' ' num2str(np(3,:))]);
    disp('Number of points in T/X/Y/Z in the plane slice:')
    disp([' ' num2str(np(4,:))]);
    disp('Number of points in T/X/Y/Z in the field beam with d/dz margins for plane-to-beam calculations:')
    disp([' ' num2str(np(5,:))]);
    disp('Number of points in T/X/Y/Z in the b2b source beam slice:')
    disp([' ' num2str(np(6,:))]);
    disp('Number of points in T/X/Y/Z in the b2b field beam slice:')
    disp([' ' num2str(np(7,:))]);
    disp('Number of points in T/X/Y/Z in the ss2b source beam slice:')
    disp([' ' num2str(np(8,:))]);
    disp('Number of points in T/X/Y/Z in the ss2b field beam slice:')
    disp([' ' num2str(np(9,:))]);
    disp('Number of points in T/X/Y/Z in the p2b source beam slice:')
    disp([' ' num2str(np(10,:))]);
    disp('Number of points in T/X/Y/Z in the p2b field beam slice:')
    disp([' ' num2str(np(11,:))]);
end

%total memory taken on each processor by the different arrays, in Megabytes
fielddistr12size = (np(1,1)+1)*np(1,2)*np(1,3)*np(1,4)*2*8/1024^2/numprocs;
sourcedistr12size = np(2,1)*np(2,2)*np(2,3)*1*2*8/1024^2/numprocs;
fieldwithmargindistr12size = (np(3,1)+1)*np(3,2)*np(3,3)*np(3,4)*2*8/1024^2/numprocs;
planedistr12size = np(4,1)*np(4,2)*np(4,3)*1*2*8/1024^2/numprocs;
fieldwithP2Bmargindistr12size = (np(5,1)+1)*np(5,2)*np(5,3)*np(5,4)*2*8/1024^2/numprocs;
b2bsourceslicedistr2size = np(6,2)*np(6,3)*np(6,4)*2*8/1024^2;
b2bfieldslicedistr2size = np(7,2)*np(7,3)*np(7,4)*2*8/1024^2;
ss2bsourceslicedistr2size = np(8,2)*np(8,3)*np(8,4)*2*8/1024^2;
ss2bfieldslicedistr2size = np(9,2)*np(9,3)*np(9,4)*2*8/1024^2;
p2bsourceslicedistr2size = np(10,2)*np(10,3)*np(10,4)*2*8/1024^2;
p2bfieldslicedistr2size = np(11,2)*np(11,3)*np(11,4)*2*8/1024^2;

%maximum used memory in Megabytes on specific instances in the program
% 1st: in the source-slice-to-beam solution, during the convolution
% 2nd: in the source-slice-to-beam solution, during the 2->1 redistribution afterwards
% 3rd: in the beam-to-beam solution, during the convolution
% 4th: in the beam-to-beam solution, during the 1->2 or 2->1 redistribution
% 5th: in the field-plane-to-beam solution, during the convolution
% 6th: in the field-plane-to-beam solution, during the 2->1 redistribution afterwards

memsize = [ ...
    sourcedistr12size + fieldwithmargindistr12size + ss2bsourceslicedistr2size + ss2bfieldslicedistr2size...
    2*fieldwithmargindistr12size ...
    fielddistr12size + b2bsourceslicedistr2size + b2bfieldslicedistr2size ...
    2*fielddistr12size...
    planedistr12size + fieldwithP2Bmargindistr12size + p2bsourceslicedistr2size + p2bfieldslicedistr2size...
    2*fieldwithP2Bmargindistr12size...
    ];

numpoints = np;

end
%% end check_ICS_memory_nomex

%% get_ICS_arraysizes

function numpoints = get_ICS_arraysizes(mediumparams,domainparams,...
    modelparams,sourceparams,numprocs,dTdelay,iDimSource)

% numpoints = get_ICS_arraysizes(mediumparams,domainparams,...
%     modelparams,sourceparams,numprocs,dTdelay,iDimSource)
% Calculates the arraysizes used internally in the ICS program

% 0- Check freq0 in different structures and correct for st -- operations
% as performed in generate_ICS_configurationfile
% correct these values directly in the different structures
freq0 = modelparams.freq0;
[lt,st] = normalize_lt_st_to_freq0(domainparams,freq0);

% - account for st in Tpulse, Tdelay and in lt
%   the employment of st has not been implemented in the Fortran program,
%   but accounting for it in lt, Tpulse and Tdelay provides an easy
%   bypass...
if(strcmp(sourceparams.definitiontype,'parametric')...
        &&~strcmp(sourceparams.signature.signaturetype,'file'))
        sourceparams.signature.Tpulse = sourceparams.signature.Tpulse + st;
        sourceparams.signature.Tdelay = sourceparams.signature.Tdelay + st;
end
lt = lt + st;
domainparams.lt = lt;
domainparams.st = 0;

% 1- Equivalent of SpecialChecksConfigParameters subroutine in Fortran
special_checks_config_parameters

% 2- Equivalent of NormalizeConfigParameters subroutine in Fortran
normalize_config_parameters

% 3- Equivalent of InitSpaceFromPhysicalDims in Fortran for RefBeam
init_space_from_physical_dims

% 4- Equivalent of InitGridSpace in Fortran for RefBeam
[domainparams.dimt,domainparams.dimx,domainparams.dimy,domainparams.dimz]...
    = init_grid_space(domainparams);

%-> domainparams now contains the dimensions for the beam space

% 5- Equivalent of InitConvolutionSpaces in Fortran for a ref beam
[ssliceref,fsliceref] = init_convolution_spaces(domainparams,domainparams,0);

%-> ssliceref and fsliceref now contain the dimensions for the xyz slices
% for a normal beam-to-beam calculation

% 6- Equivalent of operations in PrimarySourcetoField in Fortran
if(~strcmp(sourceparams.aperture.aperturetype,'pointsource')&&...
        strcmp(sourceparams.excitationtype,'fsource'))
    margins(1) = ceil(abs(modelparams.FDXorder/2 * domainparams.tant) - 1e-5);
    margins(2) = ceil(abs(modelparams.FDXorder/2 * domainparams.tanx) - 1e-5);
    margins(3) = ceil(abs(modelparams.FDXorder/2 * domainparams.tany) - 1e-5);
    margins(4) = ceil(modelparams.FDXorder/2);
else
    margins = [ 0 0 0 0 ];
end

% Equivalent of SourceDimensions in Fortran
source_dimensions

sourceparams.dimt = domainparams.dimt + 2*margins(1);
sourceparams.dimz = 1;
sourceparams.tanx = domainparams.tanx;
sourceparams.tany = domainparams.tany;

spacemargin.dimt = sourceparams.dimt;
spacemargin.dimx = domainparams.dimx + 2*margins(2);
spacemargin.dimy = domainparams.dimy + 2*margins(3);
spacemargin.dimz = domainparams.dimz + 2*margins(4);
spacemargin.tanx = domainparams.tanx;
spacemargin.tany = domainparams.tany;
[spacemargin.dimt,spacemargin.dimx,spacemargin.dimy,spacemargin.dimz]...
    = init_grid_space(spacemargin);

%-> spacemargin now contains the dimensions for the beam including d/dz
%margins for the source (if sourcetype ~=fsource, then no margins are used)

[sourceparams.dimt,sourceparams.dimx,sourceparams.dimy,sourceparams.dimz]...
    = init_grid_space(sourceparams);

%-> sourceparams now contains the dimensions for the source space

sourcetemp.dimt = sourceparams.dimt;
sourcetemp.dimx = sourceparams.dimx;
sourcetemp.dimy = sourceparams.dimy;
sourcetemp.dimz = 1;
sourcetemp.tanx = sourceparams.tanx;
sourcetemp.tany = sourceparams.tany;

% 7- Equivalent of InitConvolutionSpaces in Fortran for a beam with margin
% for the primary-source-to-field calculation
[sslicemar,fslicemar] = init_convolution_spaces(sourcetemp,spacemargin,1);

%-> sslicemar and fslicemar now contain the dimensions for the xyz slice for a beam
%with d/dz margins for the primary-source-to-beam calculations

if (numel(modelparams.beamiterations)>2)
    % 8- Equivalent of operations in CFPlaneSourcetoField in Fortran
    margins(1) = ceil(abs(modelparams.FDXorder/2 * domainparams.tant) - 1e-5);
    margins(2) = ceil(abs(modelparams.FDXorder/2 * domainparams.tanx) - 1e-5);
    margins(3) = ceil(abs(modelparams.FDXorder/2 * domainparams.tany) - 1e-5);
    margins(4) = ceil(modelparams.FDXorder/2);

    spacemarginP2B.dimt = domainparams.dimt + 2*margins(1);
    spacemarginP2B.dimx = domainparams.dimx + 2*margins(2);
    spacemarginP2B.dimy = domainparams.dimy + 2*margins(3);
    spacemarginP2B.dimz = domainparams.dimz + 2*margins(4);
    spacemarginP2B.tanx = domainparams.tanx;
    spacemarginP2B.tany = domainparams.tany;

    [spacemarginP2B.dimt,spacemarginP2B.dimx,spacemarginP2B.dimy,spacemarginP2B.dimz]...
        = init_grid_space(spacemarginP2B);

    %-> spacemarginP2B now contains the dimensions for the beam including d/dz
    %margins

    plane.dimt = spacemarginP2B.dimt;
    plane.dimx = domainparams.dimx;
    plane.dimy = domainparams.dimy;
    plane.dimz = domainparams.dimz;
    plane.tanx = domainparams.tanx;
    plane.tany = domainparams.tany;
    [plane.dimt,plane.dimx,plane.dimy,plane.dimz]...
        = init_grid_space(plane);

    %-> plane now contains the dimensions for the plane space

    planetemp = plane;
    planetemp.dimz = 1;

    % 9- Equivalent of InitConvolutionSpaces in Fortran for a beam with margin
    % for the P2B calculation
    [sslicemarP2B,fslicemarP2B] = init_convolution_spaces(planetemp,spacemarginP2B,1);

    %-> sslicemarP2B and fslicemarP2B now contain the dimensions for the xyz slice for a beam
    %with d/dz margins for the plane-to-beam calculation
else
    spacemarginP2B.dimt = 0;
    spacemarginP2B.dimx = 0;
    spacemarginP2B.dimy = 0;
    spacemarginP2B.dimz = 0;
    plane = spacemarginP2B;
    sslicemarP2B = spacemarginP2B;
    fslicemarP2B = spacemarginP2B;
end

%10- Fill numpoints with the correct sizes
numpoints = zeros(11,4);
numpoints = [...
    domainparams.dimt domainparams.dimx domainparams.dimy domainparams.dimz ;...
    sourceparams.dimt sourceparams.dimx sourceparams.dimy sourceparams.dimz ;...
    spacemargin.dimt  spacemargin.dimx  spacemargin.dimy  spacemargin.dimz  ;...
    plane.dimt        plane.dimx        plane.dimy        plane.dimz        ;...
    spacemarginP2B.dimt spacemarginP2B.dimx spacemarginP2B.dimy spacemarginP2B.dimz ;...
    ssliceref.dimt    ssliceref.dimx    ssliceref.dimy    ssliceref.dimz    ;...
    fsliceref.dimt    fsliceref.dimx    fsliceref.dimy    fsliceref.dimz    ;...
    sslicemar.dimt    sslicemar.dimx    sslicemar.dimy    sslicemar.dimz    ;...
    fslicemar.dimt    fslicemar.dimx    fslicemar.dimy    fslicemar.dimz    ;...
    sslicemarP2B.dimt sslicemarP2B.dimx sslicemarP2B.dimy sslicemarP2B.dimz ;...
    fslicemarP2B.dimt fslicemarP2B.dimx fslicemarP2B.dimy fslicemarP2B.dimz ;...
    ];

    function special_checks_config_parameters

        % if planewave, set aperture to point source
        if (strcmp(sourceparams.excitationtype,'planewave'))
            sourceparams.aperture.aperturetype='pointsource';
        end
        % y symmetry may be employed??
        if ((abs(domainparams.thetay-pi/2)<1e-5)&&...
                strcmp(modelparams.contrastsourcetype,'nonlin'))
            modelparams.useysymmetry= 1;
            domainparams.ly = domainparams.ly/2;
        else
            modelparams.useysymmetry= 0;
        end
        % set start index in z
        if(abs(domainparams.sz)<1e-6)
            domainparams.sz = mediumparams.c0 / (2*freq0*modelparams.Fnyq) * 1e3;
            domainparams.lz = domainparams.lz - domainparams.sz;
        end

    end % function special_checks_config_parameters

    function normalize_config_parameters

        lambdamm = mediumparams.c0 / freq0 * 1e3;

        % normalize dimensions
        domainparams.lx = domainparams.lx / lambdamm;
        domainparams.ly = domainparams.ly / lambdamm;
        domainparams.lz = domainparams.lz / lambdamm;
        domainparams.sz = domainparams.sz / lambdamm;
        domainparams.lt = domainparams.lt;

        % normalize source parameters
        if strcmp(sourceparams.definitiontype,'parametric')
            if strcmp(sourceparams.aperture.aperturetype,'cylindrical')
                sourceparams.aperture.radius = sourceparams.aperture.radius / lambdamm;
                sourceparams.aperture.radfocus = sourceparams.aperture.radfocus / lambdamm;
                if(abs(sourceparams.aperture.radfocus)>1e-5)
                    sourceparams.aperture.maxphasedelay = ...
                        sqrt(sourceparams.aperture.radius^2 + ...
                        sourceparams.aperture.radfocus^2)-abs(sourceparams.aperture.radfocus);
                else
                    sourceparams.aperture.maxphasedelay = 0;
                end
            elseif strcmp(sourceparams.aperture.aperturetype,'rectangular')
                sourceparams.aperture.rectwidth=sourceparams.aperture.rectwidth / lambdamm;
                sourceparams.aperture.rectheight=sourceparams.aperture.rectheight / lambdamm;
                sourceparams.aperture.maxphasedelay=0;
            elseif strcmp(sourceparams.aperture.aperturetype,'triangular')
                sourceparams.aperture.triangedge=sourceparams.aperture.triangedge / lambdamm;
                sourceparams.aperture.triangxcenter=sourceparams.aperture.triangxcenter / lambdamm;
                sourceparams.aperture.maxphasedelay=0;
            elseif strcmp(sourceparams.aperture.aperturetype,'triangular')
                sourceparams.aperture.maxphasedelay=0;
            elseif strcmp(sourceparams.aperture.aperturetype,'phasedarray')
                sourceparams.aperture.elwidth=sourceparams.aperture.elwidth / lambdamm;
                sourceparams.aperture.elheight=sourceparams.aperture.elheight / lambdamm;
                sourceparams.aperture.kerf=sourceparams.aperture.kerf / lambdamm;
                sourceparams.aperture.focusx=sourceparams.aperture.focusx / lambdamm;
                sourceparams.aperture.focusz=sourceparams.aperture.focusz / lambdamm;
                sourceparams.aperture.elevationfocusz=sourceparams.aperture.elevationfocusz / lambdamm;
                sourceparams.aperture.arraywidth=sourceparams.aperture.numels*...
                    (sourceparams.aperture.elwidth+sourceparams.aperture.kerf)-...
                    sourceparams.aperture.kerf;
                elx=-sourceparams.aperture.arraywidth/2+sourceparams.aperture.elwidth/2+ ...
                    (0:sourceparams.aperture.numels-1)*(sourceparams.aperture.elwidth+...
                    sourceparams.aperture.kerf);
                eltdelay=-sqrt((sourceparams.aperture.focusx-elx).^2+sourceparams.aperture.focusz^2);
                maxphasedelayx=max(eltdelay-min(eltdelay));

                if(sourceparams.aperture.elevationfocusz>1.0e-5)
                    maxphasedelayy = sqrt((sourceparams.aperture.elheight/2)^2 ...
                        + sourceparams.aperture.elevationfocusz^2) - abs(sourceparams.aperture.elevationfocusz);
                else
                    maxphasedelayy = 0
                end

                sourceparams.maxphasedelay = maxphasedelayx + maxphasedelayy;

            elseif strcmp(sourceparams.aperture.aperturetype,'matrixarray')

                disp('yet to be implemented')

            end
        end
    end %function normalize_config_parameters

    function init_space_from_physical_dims

        beamlz = domainparams.lz/numel(modelparams.beamiterations);

        %convert thetas to tans
        if(domainparams.thetax == pi/2)
            domainparams.tanx = 0;
        else
            domainparams.tanx = 1/tan(domainparams.thetax);
        end
        if(domainparams.thetay == pi/2)
            domainparams.tany = 0;
        else
            domainparams.tany = 1/tan(domainparams.thetay);
        end
        if(domainparams.thetat==pi/2)
            domainparams.tant = 0;
        else
            domainparams.tant = 1/tan(domainparams.thetat);
        end

        %calculate step sizes
        modelparams.dx = 1/(2*modelparams.Fnyq);
        modelparams.dt = 1/(2*modelparams.Fnyq);

        %calculate dims
        domainparams.dimt = ceil(domainparams.lt/modelparams.dt - 1e-5);
        domainparams.dimx = ceil(domainparams.lx/modelparams.dx - 1e-5);
        domainparams.dimy = ceil(domainparams.ly/modelparams.dx - 1e-5);
        domainparams.dimz = ceil(beamlz/modelparams.dx - 1e-5);

    end % function init_space_from_physical_dims

    function [dimtnew,dimxnew,dimynew,dimznew] = init_grid_space(strucin)

        dimt = strucin.dimt;
        dimx = strucin.dimx;
        dimy = strucin.dimy;
        dimz = strucin.dimz;

        %resize to 'nice' number of points (for numprocs)
        dimtnew = ceil((dimt+1)/numprocs)*numprocs - 1;
        
        if (Isprodsmallprimes(dimx*dimy*dimz) & (mod(dimx*dimy*dimz,numprocs)==0))

            dimxnew = dimx;
            dimynew = dimy;
            dimznew = dimz;

        else

            iPa                = 1;
            iPb                = 1;
            iPc                = 1;
            i                = numprocs;
            j                = 2; %first prime number
            k                = 0; %counter for the xyz-dimensions
            while (i ~= 1)
                % If not j | i, we should try the next number
                while (mod(i, j) ~= 0)
                    j = j + 1;
                end

                if (mod(k,3)==0)
                    iPa = iPa * j;
                elseif (mod(k,3)==1)
                    iPb = iPb * j;
                elseif (mod(k,3)==2)
                    iPc = iPc * j;
                end

                i        = floor(i / j);
                k         = k + 1;
            end

            divx = ceil(dimx/iPa);
            while (~Isprodsmallprimes(divx))
                divx=divx+1;
            end
            divy = ceil(dimy/iPb);
            while (~Isprodsmallprimes(divy))
                divy=divy+1;
            end
            divz = ceil(dimz/iPc);
            while (~Isprodsmallprimes(divz))
                divz=divz+1;
            end

            dimxnew = divx * iPa;
            dimynew = divy * iPb;
            dimznew = divz * iPc;

        end

    end % function init_grid_space

    function source_dimensions

        if strcmp(sourceparams.definitiontype,'file')
            sourceparams.dimt = iDimSource(1) + ceil(dTdelay/modelparams.dt);
            sourceparams.dimx = iDimSource(2);
            sourceparams.dimy = iDimSource(3);
            sourceparams.dimz = 1;
        else
            if strcmp(sourceparams.signature.signaturetype,'gaussian')
                sourceparams.dimt = ceil((sourceparams.signature.Tpulse + ...
                    sourceparams.maxphasedelay)/modelparams.dt);
            elseif strcmp(sourceparams.signature.signaturetype,'blackman')
                sourceparams.dimt = ceil((1.1 * sourceparams.signature.Tpulse + ...
                    sourceparams.maxphasedelay)/modelparams.dt);
            elseif strcmp(sourceparams.signature.signaturetype,'file')
                sourceparams.dimt = iDimSource(1) + ...
                    ceil(sourceparams.maxphasedelay/modelparams.dt);
            end
            if strcmp(sourceparams.aperture.aperturetype,'cylindrical')
                sourceparams.dimx = floor(2*(sourceparams.aperture.radius/...
                    modelparams.dx + 10) + 1);
                sourceparams.dimy = floor(2*(sourceparams.aperture.radius/...
                    modelparams.dx + 10) + 1);
            elseif strcmp(sourceparams.aperture.aperturetype,'rectangular')
                sourceparams.dimx = floor(2*(sourceparams.aperture.rectwidth/...
                    (2*modelparams.dx) + 10) + 1);
                sourceparams.dimy = floor(2*(sourceparams.aperture.rectheight/...
                    (2*modelparams.dx) + 10) + 1);
            elseif strcmp(sourceparams.aperture.aperturetype,'triangular')
                sourceparams.dimx = floor(2*(sqrt(3)*sourceparams.aperture.triangedge/...
                    (4*modelparams.dx) + abs(sourceparams.aperture.triangxcenter)/modelparams.dx + 10) + 1);
                sourceparams.dimy = floor(2*(sourceparams.aperture.triangedge/...
                    (2*modelparams.dx) + 10) + 1);
            elseif strcmp(sourceparams.aperture.aperturetype,'pointsource')
                sourceparams.dimx = 1;
                sourceparams.dimy = 1;
            elseif strcmp(sourceparams.aperture.aperturetype,'phasedarray')
                sourceparams.dimx = floor(2*(ceil(sourceparams.aperture.arraywidth/...
                    (2*modelparams.dx))+15)+1);
                sourceparams.dimy = floor(2*(ceil(sourceparams.aperture.elheight/...
                    (2*modelparams.dx))+15)+1);
            elseif strcmp(sourceparams.aperture.aperturetype,'matrixarray')
                disp('yet to be implemented');
            end
        end

    end % function source_dimensions

    function [sslice,fslice] = init_convolution_spaces(srcspace,destspace,source_is_slice)

        dimz = srcspace.dimz + destspace.dimz - 1;
        while (~Isprodsmallprimes(dimz))
            dimz = dimz + 1;
        end

        if (abs(srcspace.tanx)<1E-5 && abs(destspace.tanx)<1E-5)
            dimx = srcspace.dimx + destspace.dimx - 1;
        else
            skewxl = ceil(max([0 (srcspace.dimz * (srcspace.tanx-destspace.tanx) - 1E-5) ]));
            skewxr = ceil(max([0 (srcspace.dimz * (destspace.tanx-srcspace.tanx) - destspace.dimx - 1E-5) ]));
            dimx = srcspace.dimx + destspace.dimx - 1 + skewxl + skewxr + 2;
        end
        while (~Isprodsmallprimes(dimx))
            dimx = dimx + 1;
        end

        if (abs(srcspace.tany)<1E-5 && abs(destspace.tany)<1E-5)
            dimy = srcspace.dimy + destspace.dimy - 1;
        else
            skewyl = ceil(max([0 (srcspace.dimz * (srcspace.tany-destspace.tany) - 1E-5) ]));
            skewyr = ceil(max([0 (srcspace.dimz * (destspace.tany-srcspace.tany) - destspace.dimy - 1E-5) ]));
            dimy = srcspace.dimy + destspace.dimy - 1 + skewyl + skewyr + 2;
        end
        while (~Isprodsmallprimes(dimy))
            dimy = dimy + 1;
        end

        sslice.dimt = numprocs-1;
        sslice.dimx = dimx;
        sslice.dimy = dimy;
        if(source_is_slice==1)
            sslice.dimz = 1;
        else
            sslice.dimz = dimz;
        end

        fslice.dimt = numprocs-1;
        fslice.dimx = dimx;
        fslice.dimy = dimy;
        fslice.dimz = dimz;

    end % function init_convolution_spaces

    function Ispr=Isprodsmallprimes(Number)

        primes = [2 3 5 7];

        if (Number==0)
            Ispr = 0;
            return
        end

        for loop=1:length(primes)
            while (mod(Number,primes(loop))==0)
                Number=Number/primes(loop);
            end
        end

        if (Number==1)
            Ispr = 1;
        else
            Ispr = 0;
        end

    end % function Isprodsmallprimes

end % function get_ICS_arraysizes
