function signature = set_ICS_signature_data(signaturedata,signaturefilename,filetype)

% signature = set_ICS_signature_data(signaturedata,signaturefilename,filetype)
% Generates the source signature file and the matching source signature
% structure for a vector of signature data. 
% The signature structure contains the signature settings 
% used in the set_ICS_source_parametric function. 
% Inputs:
% signaturedata:     vector in the t dimension containing the source 
%                    signature. The discretization is assumed to be 
%                    the same as in the model parameter, and the data 
%                    is assumed to be filtered.
% signaturefilename: (optional, default = 'source_signature') filename 
%                    for the source signature file, without the extension
% filetype:          (optional, default = 'dat') File type, choose between
%                    'bin' for binary ieee-be.l64 type (check 'help fopen')
%                    or 'dat' for 16-digit ascii type (check 'help save').

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(ndims(signaturedata)~=2)
    error('signaturedata should be a vector')
end
if(nargin<2 || isempty(signaturefilename))
    signaturefilename = 'source_signature';
end
if(nargin<3 || isempty(filetype))
    filetype = 'dat';
end
if(nargin>3)
    error('Too many arguments')
end

if(strcmp(filetype,'dat'))

    % Save as an ASCII textfile, with 16-digit number size.
    % In Fortran, use format specifier "(XES25.16E3)", where 
    % X is the number of elements to be read
    %
    % File format:
    %   size(signaturedata), 2 real numbers
    %   signaturedata, numel(signaturedata) real numbers

    signaturefilename = strcat(signaturefilename,'.dat');

    datasize=size(signaturedata);

    signaturedata=signaturedata(:);

    save(signaturefilename,'-ASCII','-DOUBLE','datasize','signaturedata')

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
    %   size(signaturedata),   2 int*8 numbers
    %   signaturedata, numel(signaturedata) real*8 numbers

    signaturefilename = strcat(signaturefilename,'.bin');

    fid=fopen(signaturefilename,'w','s');

    fwrite(fid,8*length(size(signaturedata)),'integer*4');
    fwrite(fid,size(excitationdata),'integer*8');
    fwrite(fid,8*length(size(signaturedata)),'integer*4');

    fwrite(fid,8*prod(size(signaturedata)),'integer*4');
    fwrite(fid,signaturedata,'real*8');
    fwrite(fid,8*prod(size(signaturedata)),'integer*4');

    fclose(fid);

end

signature = struct('signaturetype','file',...
    'signaturefilename',signaturefilename);