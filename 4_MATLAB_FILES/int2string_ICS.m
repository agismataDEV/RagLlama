function string = int2string_ICS(integer)

% string = int2string_ICS(integer)
% Converts an integer number to a string as employed in the filenames of
% the ICS output files. The format is:
% - First character is '_' or '-' for either positive or negative integers
% - Next three characters are the integer value from 0 up to 999. Always
%   three characters, 0's are used to fill it if necessary.

% JH, version 090323

if integer<0
    plusminus = '-';
else
    plusminus = '_';
end

absinteger = abs(integer);

i1 = floor(absinteger/1000);
i2 = floor((absinteger - i1*1000)/100);
i3 = floor((absinteger - i1*1000 - i2*100)/10);
i4 = floor(absinteger - i1*1000 - i2*100 - i3*10);

string = [ plusminus int2str(i2) int2str(i3) int2str(i4) ]; 
