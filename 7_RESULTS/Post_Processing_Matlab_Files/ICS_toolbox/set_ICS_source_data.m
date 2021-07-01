function sourceparams = set_ICS_source_data(excitationdata,excitationtype,...
    filename,filetype,timedelaydata)

% sourceparams = set_ICS_source_data(excitationdata,excitationtype,...
%     filename,filetype,timedelaydata)
% Generates the source data file and the matching sourceparams structure
% containing all source settings used in the ICS Parnac program.
% Inputs:
% excitationdata: three-dimensional array in the t,x and y dimensions
%                 containing the source plane excitation. It is assumed 
%                 that the step sizes in t, x and y are correct, that 
%                 the data is centered around (x,y) = (0,0), and that the 
%                 data is assumed to be filtered.
% excitationtype: (optional, default = 'fsource') Type of source excitation,
%                 choose between 'qsource','fsource' or 'dtqsource'
% filename:       (optional, default = 'source_excitation') filename 
%                 for the source data file, without the extension
% filetype:       (optional, default = 'bin') File type, choose between
%                 'bin' for binary ieee-be.l64 type (check 'help fopen')
%                 or 'dat' for 16-digit ascii type (check 'help save').
%                 WARNING: timedelaydata is not stored in 'dat' data type!!
% timedelaydata:  (optional, default is zero-filled array) two-dimensional
%                 array in the x and y dimensions containing the time
%                 delays in [ms].

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(ndims(excitationdata)~=3)
    error('excitationdata should have three dimensions')
end
if(nargin<2 || isempty(excitationtype))
    excitationtype = 'fsource';
end
if(nargin<3 || isempty(filename))
    filename = 'source_excitation';
end
if(nargin<4 || isempty(filetype))
    filetype = 'bin';
end
if(nargin<5 || isempty(timedelaydata))
    timedelaydata=zeros(size(excitationdata,2),size(excitationdata,3));
end
if(nargin>5)
    error('Too many arguments')
end

if(isempty(strmatch(excitationtype,strvcat(...
        'qsource','fsource','dtqsource'), 'exact')))
    error('excitationtype has an invalid value')
end

if(isempty(strmatch(filetype,strvcat(...
        'bin','dat'), 'exact')))
    error('filetype has an invalid value')
end

if( (size(timedelaydata,1)~=size(excitationdata,2))...
        || (size(timedelaydata,2)~=size(excitationdata,3)))
    error('Size of timedelaydata does not match size of excitationdata')
end

if(strcmp(filetype,'dat'))

    % Save as an ASCII textfile, with 16-digit number size.
    % In Fortran, use format specifier "(XES25.16E3)", where 
    % X is the number of elements to be read
    %
    % File format:
    %   size(excitationdata), 3 real numbers
    %   excitationdata, numel(excitationdata) real numbers, sorted in
    %       column-major order

    filename = strcat(filename,'.dat');

    datasize=size(excitationdata);

    excitationdata=excitationdata(:);

    save(filename,'-ASCII','-DOUBLE','datasize','excitationdata')

elseif(strcmp(filetype,'bin'))

    % Save as a Fortran unformatted binary file.
    % data type: IEEE floating point with big-endian byte ordering
    %               and 64 bit long.
    %
    % Beware: in the Fortran unformatted binary file,
    % before and after each entry the number of bytes used by the entry
    % must be given as an integer*4!!
    %
    % File format:
    %   size(excitationdata),   3 int*8 numbers
    %   max(max(timedelaydata), 1 real*8 number
    %   timedelaydata, numel(timedelaydata) real*8 numbers,
    %       sorted in column-major order
    %   excitationdata, numel(excitationdata) real*8 numbers,
    %       sorted in column-major order
    

    filename = strcat(filename,'.bin');

    fid=fopen(filename,'w','s');

    fwrite(fid,8*length(size(excitationdata)),'integer*4');
    fwrite(fid,size(excitationdata),'integer*8');
    fwrite(fid,8*length(size(excitationdata)),'integer*4');

    fwrite(fid,8,'integer*4');
    fwrite(fid,max(max(timedelaydata)),'real*8');
    fwrite(fid,8,'integer*4');

    fwrite(fid,8*prod(size(timedelaydata)),'integer*4');
    fwrite(fid,timedelaydata,'real*8');
    fwrite(fid,8*prod(size(timedelaydata)),'integer*4');

    fwrite(fid,8*prod(size(excitationdata)),'integer*4');
    fwrite(fid,excitationdata,'real*8');
    fwrite(fid,8*prod(size(excitationdata)),'integer*4');

    fclose(fid);

end

sourceparams = struct('excitationtype',excitationtype,...
    'definitiontype','file',...    
    'filename',filename);
