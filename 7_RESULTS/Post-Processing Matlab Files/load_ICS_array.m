function output = load_ICS_array(filename)

% output = load_ICS_array(filename)
% Reads arrays as exported in HDF5 format in the ICS Parnac program 
% by the procedure ExportArray.
% Input:
% filename: full name of the file, including extension
% Output:   the structure output with the fields
% data:     containing the data, either real or complex 
%           and with 1...3 dimensions.

% JH - version 090126

fileinfo = hdf5info(filename);
dataset = fileinfo.GroupHierarchy.Datasets(1);

% Read data
output.data=hdf5read(filename,dataset.Name);

% If the stored array is complex, read the second stored array 
% as the imaginary term
if (length(fileinfo.GroupHierarchy.Datasets)==2)
    dataset = fileinfo.GroupHierarchy.Datasets(2);
    output.data=output.data+1i*hdf5read(filename,dataset.Name);
end

% The first two dimensions of the array seem to be permuted. Permute them
% back...
if (ndims(output.data)==3)
    output.data=permute(output.data,[2 1 3]);
elseif (ndims(output.data)==2)
    output.data=permute(output.data,[2 1]);
else
    disp('load_ICS_array: dimensions may be permuted...');
end
