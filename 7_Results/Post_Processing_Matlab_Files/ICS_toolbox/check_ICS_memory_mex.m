function [memsize,numpoints] = check_ICS_memory_mex(configfilename,numprocs,showoutput,dTdelay,iDimSource)

% [memsize,numpoints] = check_ICS_memory_mex(configfilename,numprocs,showoutput,dTdelay,iDimSource)
% Calculates the memorysizes per processor for several critical instances in the 
% ICS program. Whereas check_ICS_memory works with a Matlab emulation of the Fortran code,
% check_ICS_memory_mex works with a Matlab executable (mex) that
% directly evaluates parts of the Fortran program.
% Inputs:
% configfilename: name of configuration file, to be found in the current directory
% numprocs:       number of processors
% showoutput:     (optional, default = 0) echo the configuration file (1) or not (0)
% dTdelay:        (optional, default = 0) dTdelay of the source data input file
%                 if any is used
% iDimSource:     (optional, default = [0 0 0]) vector of three numbers giving 
%                 the T/X/Y dimension sizes of the source data input file. 
%                 If a source signature input file is used, only the T 
%                 dimension is used, the rest is derived from the parametric 
%                 aperture.
% Outputs:
% memsize:        vector with the critical memorysizes per processor in MB
% numpoints:     11x4 element array containing the number of points in the
%                T/X/Y/Z dimensions. The rows correspond with the
%                following arrays employed internally in the ICS program:
%                1 - The field beam
%                2 - The source slice
%                3 - The field  beam including d/dz margins for source2beam
%                4 - The plane slice
%                5 - The field  beam including d/dz margins for plane2beam
%                6 - The source beam slice for the beam-to-beam calculations
%                7 - The field  beam slice for the beam-to-beam calculations
%                8 - The source beam slice for the source-slice-to-beam calculations
%                9 - The field  beam slice for the source-slice-to-beam calculations
%                10 - The source beam slice for the plane-to-beam calculations
%                11 - The field  beam slice for the plane-to-beam calculations
%
% The function needs the mex-file ParnacRunEstimator.mexw32 to be located
% somewhere in the path.

if(nargin<2)
    error('Not enough arguments')
end
if(nargin<3 || isempty(showoutput))
    showoutput = 0;
end
if(nargin<4 || isempty(dTdelay))
    dTdelay = 0;
end
if(nargin<5 || isempty(iDimSource))
    iDimSource = [0 0 0];
end
if(nargin>5)
    error('Too many arguments')
end

np = ParnacRunEstimator(configfilename,numprocs,showoutput,dTdelay,iDimSource);

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
    2*fieldwithP2Bmargindistr12size
    ];

numpoints = np;
