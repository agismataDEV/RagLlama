function generate_ICS_LLscript(filename,numprocs,configurationfilename,...
    outputrootname,inputdir,outputdir,timeestimate)

% generate_ICS_LLscript(filename,numprocs,configurationfilename,
%       outputrootname,inputdirname,outputdirname,timeestimate)
% Generates the run script for multi-processor computer employing
% the LoadLeveler job management program.
% Inputs:
% filename:              filename of the script
% numprocs:              number of processors
% configurationfilename: filename of the configuration file
% outputrootname:        (optional, default = root of configuration file)
%                        root name of the output files
% inputdir:              (optional, default = /home/parnac/input/) path of the input directory
% outputdir:             (optional, default = /home/parnac/output/) path of the output directory        
% timeestimate:          (optional, default = 40) runtime estimate
%                        in [hr]

% JH, version 090116

if(nargin<3)
    error('Not enough arguments')
end
if(nargin<4||isempty(outputrootname))
    k=strfind(configurationfilename,'.');
    if(isempty(k))
        outputrootname = configurationfilename;
    else
        outputrootname = configurationfilename(1:k(numel(k))-1);
    end
end
if(nargin<5||isempty(inputdir))
    inputdir = '/home/parnac/input/';
end
if(nargin<6||isempty(outputdir))
    outputdir = '/home/parnac/output/';
end
if(nargin<7||isempty(timeestimate))
    timeestimate = 40;
end
if(nargin>7)
    error('Too many arguments')
end

% test numproc: <32, arbitrary, or >32 and then x*32
if((numprocs>31)&&(mod(numprocs,32)~=0))
    error('If numproc>32, then numproc should be a multiple of 32')
end

%translate runtime
hours = floor(timeestimate);
mins = floor((timeestimate-hours)*60);
secs = floor((timeestimate - hours)*3600 - mins*60);
if(hours>99)
    error('runtime should be below 100 hours')
end

fid=fopen(filename,'w');

% Remark: Fortran uses only \n as end-of-line character
% Windows would need \r\n as end-of-line

if(numprocs<32)
    fprintf(fid,'# This job will run on 1 shared node and use X processors, X<32\n');
    fprintf(fid,'# @ node = 1\n');
    fprintf(fid,'# \n');
    fprintf(fid,'# request X processors, and X hours\n');
    fprintf(fid,'#\n');
    fprintf(fid,'# @ tasks_per_node = %2.0f\n',numprocs);
else
    fprintf(fid,'# This job will run on X nodes, 32 processors on each node,\n');
    fprintf(fid,'# making a total of X*32 processors.\n');
    fprintf(fid,'# @ node = %2.0f\n',floor(numprocs/32));
    fprintf(fid,'# @ tasks_per_node = 32\n');
end

fprintf(fid,'# Job will take at most X hours wallclock time\n');
fprintf(fid,'# @ wall_clock_limit = %02.0f:%02.0f:%02.0f\n',hours,mins,secs);
fprintf(fid,'#\n');
fprintf(fid,'# standard settings:\n');
fprintf(fid,'#\n');
fprintf(fid,'# never send notification email\n');
fprintf(fid,'# @ notification = never\n');
fprintf(fid,'# define stdin, stdout and stderr\n');
fprintf(fid,'# @ input = /dev/null\n');
fprintf(fid,'# @ output = run.$(jobid).out\n');
fprintf(fid,'# @ error = run.$(jobid).err\n');
fprintf(fid,'# Tell Loadleveler that this is a parallel job\n');
fprintf(fid,'# @ job_type = parallel\n');
if(numprocs<32)
    fprintf(fid,'# This job will run on a shared node\n');
    fprintf(fid,'# @ node_usage = shared			\n');
else
    fprintf(fid,'# communication between nodes use the infiniband hardware\n');
    fprintf(fid,'# @ network.MPI = sn_all,not_shared,US	\n');
end
fprintf(fid,'# @ queue\n');
fprintf(fid,'#\n');
fprintf(fid,'# Here the shell script starts. \n');
fprintf(fid,'cd $HOME/Bin\n');
fprintf(fid,['./ParnacMain ' configurationfilename ' ' outputrootname ' ' inputdir ' ' outputdir '\n']);

fclose(fid);
