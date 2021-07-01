function struc = load_parnac_slice(filename)

% This function can be used to load hdf5 output files
% from the parnac program - JH, version 070407

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
    return
end

%disp('Importing data...')
drawnow

fileinfo = hdf5info(totalfilename);
dataset = fileinfo.GroupHierarchy.Datasets(1);

[data,attributes]=hdf5read(totalfilename,dataset.Name,'readAttributes', true);

if(max(isnan(data))>.5)
    disp(['Data contains NaN in ' totalfilename])
end

attributes_cell = {attributes.Shortname};
%%
iStart=double(attributes(strcmp(attributes_cell, 'iStart')).Value)';
iToffset=double(attributes(strcmp(attributes_cell, 'iToffset')).Value)';
iXoffset=double(attributes(strcmp(attributes_cell, 'iXoffset')).Value)';
iYoffset=double(attributes(strcmp(attributes_cell, 'iYoffset')).Value)';
rTheta=double(attributes(strcmp(attributes_cell, 'rTheta')).Value)';
rStepsize=double(attributes(strcmp(attributes_cell, 'rStepsize')).Value)';
iSaveLength=double(attributes(strcmp(attributes_cell, 'iSaveLength')).Value)';
iPermute=double(attributes(strcmp(attributes_cell, 'iPermute')).Value)';
iNumProc=double(attributes(strcmp(attributes_cell, 'iNumProc')).Value)';
iIsReal=double(attributes(strcmp(attributes_cell, 'iIsReal')).Value)';
 
%se fosse un dato complesso avrei bisogno del doppio della dimensione epr
%il

if (iIsReal)
    datatot=zeros(prod(iSaveLength),1);
else
    datatot=zeros(2*prod(iSaveLength),1);
end
dataind=length(data);
datatot(1:dataind)=data';

%firstproc=0;

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

if (~iIsReal)
    datatot=datatot(1:2:length(datatot))+1i*datatot(2:2:length(datatot));
end

%disp('Reshaping data...')
drawnow
% keyboard
datatot=reshape(datatot,iSaveLength);
datatot=permute(datatot,iPermute);

struc.start=iStart;
struc.TXYoff=[iToffset' iXoffset' iYoffset'];
struc.stepsize=rStepsize;
struc.data=datatot;
