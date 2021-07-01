% make_INCS_1MHz_48el_steered

% make the configuration file and LoadLeveler runscript for a 48 element phased
% array transducer transmitting at 1MHz and steered at 45 degrees

% set configuration: medium, domain and model
mediumparams = set_ICS_medium(998,1492,3.52);
domainparams = set_ICS_domain(35,30,16,60,0,0,0.61547970867039,[pi/4 pi/2],1.0e6);
modelparams = set_ICS_model(1.0e6,2.5,[3 3 3],'nonlin','firstandlast','y',0);

% set parametric source aperture and signature
sourceaperture = set_ICS_aperture_phasedarray(48,.21,12,.29,40,40,40);
sourcesignature = set_ICS_signature_gaussian(2.5e5,12,3,6,2,1.0e6);
sourceparams = set_ICS_source_parametric(sourceaperture,sourcesignature, 'fsource');

% set runtime environment and file/directory names
numprocs = 4;
configfilename = 'config_INCS_1MHz_48el_steered.in';
runscriptname = 'run_INCS_1MHz_48el_steered';
outputroot = 'output_INCS_1MHz_48el_steered';
timeestimate = 2;

% generate configuration file and Huygens script
generate_ICS_configurationfile(configfilename,mediumparams,domainparams,...
    modelparams,sourceparams)

generate_ICS_LLscript(runscriptname,numprocs,configfilename,...
    outputroot,[],[],timeestimate)

% check size of the timeframe
[lt_opt,thetat_opt,st_opt] = check_ICS_timeframe(mediumparams,domainparams,...
    modelparams,sourceparams,1);

% check memory usage
[memsize,numpoints] = check_ICS_memory(mediumparams,domainparams,...
    modelparams,sourceparams,numprocs,0);

disp(sprintf(' Maximum memory per processor in MB: %0.8g',max(memsize)));
if (max(memsize)>3*1024)
    disp(' Warning: maximum memory per processor is larger than 3GB; this may give problems')
end

% estimate wall clock time
wall_clock_time = check_ICS_runtime(modelparams,numpoints, numprocs, mediumparams.a_att);
disp(sprintf(' Estimated wall clock time in hours: %0.8g',wall_clock_time));
disp(sprintf(' Estimated total node hours: %0.8g',wall_clock_time*numprocs));
if(wall_clock_time>timeestimate)
    disp(' Warning: estimated wall clock time is larger than the timeestimate setting')
end

