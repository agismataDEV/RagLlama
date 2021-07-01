% test_ICS_toolbox

% test the different functions in the ICS toolbox

% set configuration: medium, domain and model
mediumparams = set_ICS_medium(1.0e3,1.0e3,3.5);
domainparams = set_ICS_domain(50,23,38,60,0,0,pi/4,[pi/2 pi/8],1.0e6);
modelparams = set_ICS_model(1.0e6,2.5,[1],'nonlin','firstandlast',...
    'y',0,4,[true true true],[30 12]);

% set parametric source aperture and signature
sourceaperture = set_ICS_aperture_phasedarray(64,.25,16,.05,-20,40,60);
% sourceaperture = set_ICS_aperture_phasedarray(64,.25,16,.05,40,40,60,...
%   .5*ones(1,64),'apodi');
% sourceaperture = set_ICS_aperture_rectangular(20,15);
% sourceaperture = set_ICS_aperture_cylindrical(20,10);
% sourceaperture = set_ICS_aperture_pointsource;
sourcesignature = set_ICS_signature_gaussian(5e5,12,3,6,4,1.0e6);
sourceparams = set_ICS_source_parametric(sourceaperture,sourcesignature);

% set source data file
% excitationdata = zeros(3,3,3);
% excitationdata(1,:,:)=1;
% excitationdata(2,:,:)=2;
% excitationdata(3,:,:)=3;
% timedelaydata = .1*ones(3,3,3);
% sourceparams = set_ICS_source_data(excitationdata,'fsource','test_source','dat',timedelaydata);

% set runtime environment and file/directory names
numprocs = 4;
configfilename = 'config_test.in';
% runscriptname = 'runscript_test2';
% outputroot = 'test';
% timeestimate = 1;

% generate configuration file and Huygens script
generate_ICS_configurationfile(configfilename,mediumparams,domainparams,...
    modelparams,sourceparams)
% 
% generate_ICS_huygensscript(runscriptname,numprocs,configfilename,...
%     outputroot,[],[],timeestimate,copyon)

% % check memory usage and timeframe
% check_ICS_memory(configfilename,numprocs,1,0,size(excitationdata))

check_ICS_memory(configfilename,numprocs,1)

 check_ICS_memory_nomex(mediumparams,domainparams,...
     modelparams,sourceparams,numprocs,1)

% [lt_opt,thetat_opt,st_opt]=check_ICS_timeframe(mediumparams,domainparams,...
%     modelparams,sourceparams,1)