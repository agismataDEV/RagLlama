function output = load_ICS_slice(filename,strides,offsets)

% output = load_ICS_slice(filename,strides,offsets)
% Load the slices as output in HDF5 format in the ICS Parnac program
% by the procedure ExportSlice.
% Inputs:
% filename: root of the filename without extension and without the
%           processor counter in the filename.
% strides:  (optional, default = [1 1 1 1]) strides that are employed in
%           each dimension to store only a limited part of the slice into
%           the output structure.
% offsets:  (optional, default = [0 0 0 0]) offsets in each dimension in
%           the global grid.
% With strides and offsets, the points (offsets + [k l m n].*strides)
% are selected out of the global grid, where [k l m n] are indices in t, x, y and z.
% Outputs:  the structure output containing the following fields:
% start:    index in the global grid of the starting point in each 
%           dimension T/X/Y/Z; multiply with the stepsize to get the
%           starting position
% TXYoff:   array of (num Z points) x 3 points that give the offset in
%           points for in the T/X/Y dimensions as a function of the 
%           Z-dimension.
% stepsize: step sizes in each of the T/X/Y/Z dimensions in [s] or [m]
% data:     The data as a 4-dimensional array in T/X/Y/Z
%

% JH - version 090126

% argument checking
if(nargin<1)
    error('Error: minimum number of arguments is 1')
end
if(nargin<2)
    strides=ones(1,4);
end
if(nargin<3)
    offsets=zeros(1,4);
end
if(nargin>3)
    error('Error: maximum number of arguments is 3')
end

if (numel(strides)~=4)|(min(strides)<1)|...
        (floor(strides)~=strides)
    error('Error: strides should be 4 integer values of 1 or larger')
end
if (numel(offsets)~=4)|(min(offsets)<0)|(max(offsets-strides)>=0)|...
        (floor(offsets)~=offsets)
    error('Error: offsets should be 4 integer values between 0 and strides-1')
end

% Locate the first file

fileexists=0;
for firstproc=0:200
    inum=['000' num2str(firstproc)];
    inum=inum(length(inum)+(-2:0));
    totalfilename=[filename '_' inum '.h5'];
    if(exist(totalfilename,'file'))
        fileexists=1;
        break
    end
end

if (fileexists==0)
    error(['cannot find files ' filename ]);
end

disp(['Importing data from ... ',filename])
drawnow

% Import the first file
fileinfo = hdf5info(totalfilename);
dataset = fileinfo.GroupHierarchy.Datasets(1);

[data,attributes]=hdf5read(totalfilename,dataset.Name,'readAttributes', true);

if(any(isnan(data)))
    disp(['Data contains NaN in ' totalfilename])
end


attributes_cell = {attributes.Shortname};
iStart=double(attributes(strcmp(attributes_cell, 'iStart')).Value)';
iToffset=double(attributes(strcmp(attributes_cell, 'iToffset')).Value)';
iXoffset=double(attributes(strcmp(attributes_cell, 'iXoffset')).Value)';
iYoffset=double(attributes(strcmp(attributes_cell, 'iYoffset')).Value)';
rTheta=double(attributes(strcmp(attributes_cell, 'rTheta')).Value)';
rStepsize=double(attributes(strcmp(attributes_cell, 'rStepsize')).Value)';
iSaveLength=double(attributes(strcmp(attributes_cell, 'iSaveLength')).Value)';
if size(iSaveLength,2) == 5 ;iSaveLength=double(attributes(strcmp(attributes_cell, 'iSaveLength')).Value(1:end-1))'; ...
        iLengthSlice=double(attributes(strcmp(attributes_cell, 'iSaveLength')).Value(end))';else iLengthSlice = 0; end
iPermute=double(attributes(strcmp(attributes_cell, 'iPermute')).Value)';
iNumProc=double(attributes(strcmp(attributes_cell, 'iNumProc')).Value)';
if sum(strcmp(attributes_cell, 'iSliceIndex')) ; iSliceIndex=double(attributes(strcmp(attributes_cell, 'iSliceIndex')).Value)'; ...
        else iSliceIndex = 0 ; end
iIsReal=double(attributes(strcmp(attributes_cell, 'iIsReal')).Value)';

iTXYoffset=[iToffset' iXoffset' iYoffset'];

if (all(strides==[1 1 1 1])&&all(offsets==[0 0 0 0]))
    %The simple case: load all data into the output structure

    if (iIsReal==1)
        datatot=zeros(prod(iSaveLength),1);
    else
        datatot=zeros(2*prod(iSaveLength),1);
    end
    dataind=length(data);
    datatot(1:dataind)=data';

    for proc=firstproc+1:iNumProc-1

        inum=['000' num2str(proc)];
        inum=inum(length(inum)+(-2:0));
        totalfilename=[filename '_' inum '.h5'];

        if (exist(totalfilename,'file'))
            fileinfo = hdf5info(totalfilename);
            dataset = fileinfo.GroupHierarchy.Datasets(1);

            data=hdf5read(totalfilename,dataset.Name);

            if(max(isnan(data))>.5)
                disp(['Data contains NaN in ' totalfilename])
            end

            datatot(dataind+(1:length(data)))=data;
            dataind=dataind+length(data);

        end
    end

    if (iIsReal~=1)
        datatot=datatot(1:2:length(datatot))+1i*datatot(2:2:length(datatot));
    end

    disp('Reshaping data...')
    drawnow

    datatot=reshape(datatot,iSaveLength);
    datatot=permute(datatot,iPermute);

    output.start=iStart;
    output.TXYoff=iTXYoffset;
    output.stepsize=rStepsize;
    output.sliceindex=iSliceIndex;
    output.lengthslice = iLengthSlice;
    output.iNumProc = iNumProc;
    output.data=datatot;

else
    %The difficult case: employing strides and offsets

    % check if distribution is 0; otherwise,
    % this function doesn't work yet...
    if ((any(iPermute~=[1 2 3 4])||iIsReal~=1))
        error('Error: load_ICS_slice with strides and offsets only works for data in Dist. 0')
    end

    % Create Corrected iSaveLength, iStart and iTXYoffset
    iPartStart = iStart;
    iPartSavelength = iSaveLength;
    addx=zeros(4,1); %the addx variable is a crucial temporary variable since
    % it decides the offset inside the beam grid in each dimension.
    % It differs from the offsets variable in that it is inside the
    % beam grid, not relative to the absolute grid as the offsets variable is.
    % Create PartStart and PartSavelength of Z-axis first
    addx(4) = mod(offsets(4) - iStart(4),strides(4));
    iPartStart(4) = iStart(4) + addx(4);
    iPartSavelength(4) = ceil((iSaveLength(4) - addx(4))/strides(4));
    % Create PartStart of TXY axes - first apply TXYoffset because of possible shift
    % in Z, then apply the same procedure as for Z
    iPartStart(1:3) = iStart(1:3) + iTXYoffset(1+addx(4),1:3);
    for TXYind=1:3
        addx(TXYind) = mod(offsets(TXYind) - iPartStart(TXYind),strides(TXYind));
        iPartStart(TXYind) = iPartStart(TXYind) + addx(TXYind);
    end
    % Create PartTXYoffset and PartSavelength of TXY axes
    % use a PartSavelength and a PartMinSavelength because the length may differ
    % on each z-location depending on iTXYoffset and number of points...
    iPartTXYoffset=zeros(iPartSavelength(4),3);
    iPartMinSavelength = iPartSavelength;
    for iPartindz=1:iPartSavelength(4)
        indz=1+addx(4)+(iPartindz-1)*strides(4);
        iPartTXYoffset(iPartindz,:) = iTXYoffset(indz,:)-(iPartStart(1:3)-iStart(1:3));
        for TXYind=1:3
            if(strides(TXYind)>1)
                addx(TXYind) = mod( offsets(TXYind) - ...
                    (iPartStart(TXYind) + iPartTXYoffset(iPartindz,TXYind))...
                    ,strides(TXYind));
                iPartMinSavelength(TXYind)=min(iPartMinSavelength(TXYind),...
                    ceil((iSaveLength(TXYind)-addx(TXYind))/strides(TXYind)));
                iPartSavelength(TXYind)=max(iPartMinSavelength(TXYind),...
                    ceil((iSaveLength(TXYind)-addx(TXYind))/strides(TXYind)));
                iPartTXYoffset(iPartindz,TXYind)=iPartTXYoffset(iPartindz,TXYind)+...
                    addx(TXYind);
            end
        end
    end

    % create data2 array
    % if (iIsReal)
    data2=zeros(prod(iPartSavelength),1);
    % else
    %     data2=zeros(2*prod(iPartSavelength),1);
    % end

    % only insert those data points of the first file into data2 which are on the
    % selected points
    xyzglobalind = 0;     %current XYZ position in data array over all files
    xyzpartind = 0; %current XYZ position in data2 array
    xyzlength=length(data)/iSaveLength(1); %total XYZ positions in this file
    for xyzlocalind=0:xyzlength-1
        xyzind = xyzglobalind+xyzlocalind;
        xind = mod(xyzind,iSaveLength(2));
        yind = mod(floor(xyzind/iSaveLength(2)),iSaveLength(3));
        zind = floor(floor(xyzind/iSaveLength(2))/iSaveLength(3));

        %include time trace only if it is on a 'good' position
        if ((mod(iStart(4) + zind + offsets(4),strides(4))==0)&&...
                (mod(iStart(3) + yind + iTXYoffset(1+zind,3) + offsets(3),strides(3))==0)&&...
                (mod(iStart(2) + xind + iTXYoffset(1+zind,2) + offsets(2),strides(2))==0))
            addx(1) = mod(offsets(1) - (iStart(1) + iTXYoffset(1+zind,1)),strides(1));
            data2(1 + iPartSavelength(1)*xyzpartind + (0:iPartSavelength(1)-1)) = ...
                data(1 + iSaveLength(1)*xyzlocalind + addx(1) + strides(1)*(0:iPartSavelength(1)-1));
            % for xyzlocalind=xyzlength-1, we may have an exceeding of the
            % data array if the last trace is iPartMinSavelength(1) instead of
            % iPartSavelength(1)... hardly ever happens, as we most often
            % take strides(1)=1!
            xyzpartind = xyzpartind+1;
        end
    end
    xyzglobalind = xyzglobalind+xyzlength;

    % read all other files and insert selected data points into data2
    for proc=firstproc+1:iNumProc-1

        inum=['000' num2str(proc)];
        inum=inum(length(inum)+(-2:0));
        totalfilename=[filename '_' inum '.h5'];

        if (exist(totalfilename,'file'))
            fileinfo = hdf5info(totalfilename);
            dataset = fileinfo.GroupHierarchy.Datasets(1);

            data=hdf5read(totalfilename,dataset.Name);

            if(any(isnan(data)))
                disp(['Data contains NaN in ' totalfilename])
            end

            %reduce the dataset in each dimension
            xyzlength=length(data)/iSaveLength(1); %total XYZ positions in this file
            for xyzlocalind=0:xyzlength-1
                xyzind = xyzglobalind+xyzlocalind;
                xind = mod(xyzind,iSaveLength(2));
                yind = mod(floor(xyzind/iSaveLength(2)),iSaveLength(3));
                zind = floor(floor(xyzind/iSaveLength(2))/iSaveLength(3));

                %include time trace only if it is on a 'good' position
                if ((mod(iStart(4) + zind + offsets(4),strides(4))==0)&&...
                        (mod(iStart(3) + yind + iTXYoffset(1+zind,3) + offsets(3),strides(3))==0)&&...
                        (mod(iStart(2) + xind + iTXYoffset(1+zind,2) + offsets(2),strides(2))==0))
                    addx(1) = mod(iStart(1) + iTXYoffset(1+zind,1) + offsets(1),strides(1));
                    data2(1 + iPartSavelength(1)*xyzpartind + (0:iPartSavelength(1)-1)) = ...
                        data(1 + iSaveLength(1)*xyzlocalind + addx(1) + strides(1)*(0:iPartSavelength(1)-1));
                    xyzpartind = xyzpartind+1;
                end
            end
            xyzglobalind = xyzglobalind+xyzlength;

        end
    end

    % if (~iIsReal)
    %     data=data2(1:2:length(data2))+1i*data2(2:2:length(data2));
    %     data2=data;
    % end

    disp('Reshaping data...')
    drawnow

    %reshape, truncate and permute data array
    data2=reshape(data2,iPartSavelength);
    data2=data2(1:iPartMinSavelength(1),1:iPartMinSavelength(2),...
        1:iPartMinSavelength(3),1:iPartMinSavelength(4));
    % data2=permute(data2,iPermute);

    %divide/multiply variables by strides and store in outputture
    for TXYind=1:4
        if(TXYind<4)
            iPartTXYoffset(:,TXYind)=iPartTXYoffset(:,TXYind)/strides(TXYind);
        end
        iPartStart(TXYind)=iPartStart(TXYind)/strides(TXYind);
        rStepsize(TXYind)=rStepsize(TXYind)*strides(TXYind);
    end
    output.start=iPartStart;
    output.TXYoff=iPartTXYoffset;
    output.stepsize=rStepsize;    
    output.sliceindex=iSliceIndex;
    output.lengthslice = iLengthSlice;
    output.data=data2;

end
