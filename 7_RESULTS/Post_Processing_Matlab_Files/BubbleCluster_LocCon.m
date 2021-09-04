% ***********************************************************************************************************
% LOADING THE LOCATION AND THE SCATTERED PRESSURE OF EACH SCATTERER
%
% **********************************************************
% [Bubble] = BubbleCluster_LocCon(medium,domain,file,Bubble)
% **********************************************************
% Plot the domain with the existence of the scatterers
% Plot all the sideviews of the domain to understand exactly where the scatterers
% are placed
%
% Inputs:
% medium        struct      Contains all the medium properties, e.g Density, Speed of Sound etc
% domain        struct      Contains all the computational domain properties, e.g space and time step
% file          struct      Contains all the loaded file properties ,e.g directory and root name etc
% Bubble        struct      Contains all properties of the Scatterers' Cluster
%
% Outputs:
% Bubble       struct      Contains all properties of the Scatterers' Cluster
%
% ******************
% AM, version 200723
% ******************
%
% ***********************************************************************************************************


function [Bubble] = BubbleCluster_LocCon(medium,domain,file,Bubble)
%%
disp(' ')
disp('Loading the microbubble cluster location and contrast ...')
if (Bubble.N >0)
    Bubble.LocFile  = zeros(2*Bubble.N,3);
    Bubble.ConFile  = zeros(2*Bubble.N,domain.tdimpar); 
    Bubble.RadFile  = zeros(2*Bubble.N,1); 
    
    % Load the bubbles' location from file
    Bubble.LocNorm  = LoadContrastParams(Bubble,file,'Bubble_Location', Bubble.LocFile , 'yes');
    % Load the bubbles' contrast source term value
    Bubble.Contrast = LoadContrastParams(Bubble,file,'Bubble_Contrast', Bubble.ConFile , file.load_contrast_from_file, 'rb');
    % Load the bubbles' radius value
    Bubble.Radius   = LoadContrastParams(Bubble,file,'Bubble_Radius'  , Bubble.RadFile , file.load_radius_from_file);
    
    % Find Location of each microbubble
    Bubble.Loc = Bubble.LocNorm*medium.dLambdaNN;
    
    Bubble.LocGlob(:,1) = Bubble.Loc(:,1) + (domain.TXYstart(2)+domain.TXYoff(1,2)) * domain.dxpar;
    Bubble.LocGlob(:,2) = Bubble.Loc(:,2) + (domain.TXYstart(3)+domain.TXYoff(1,3)) * domain.dypar;
    Bubble.LocGlob(:,3) = Bubble.Loc(:,3) + domain.start(3)*domain.dzpar;
    % Find the range of the microbubble cluster , calculated with the
    % dimensions specified in the INCS code
    Bubble.LocRange = [ min(Bubble.LocGlob); max(Bubble.LocGlob)]';
    Bubble.Volume = prod(max(Bubble.LocGlob) - min(Bubble.LocGlob));
    %% CHECK IF ANY BUBBLE IS LOCATED TOO CLOSE WITH ANOTHER BUBBLE
%     if (Bubble.N>1)
%         disp('Calculating minimum distance between Bubbles ...')
%         q = ones(Bubble.N,1);
%         for iBubble = 2:Bubble.N-1
%             q(1:iBubble-1) = sqrt(sum((Bubble.LocGlob(iBubble,:)- Bubble.LocGlob(1:iBubble-1,:)).^2,2));
%             
%             %         if (~sum(q==0)) ;break; end
%             Bubble.diff(iBubble-1) =  min(q(1:iBubble-1)) ;
%             if (Bubble.diff(iBubble-1)*1e-3< Bubble.mindist) ;error(['Error , Bubble No ' , num2str(iBubble),' is at ', num2str(min(Bubble.diff*1e-3)), ...
%                     '[m] distance to a neighboring bubble']) ; end
%             
%         end
%         [MinVal , MinIndex] = min(Bubble.diff(1:end-1)*1e-3);
        
%         q = ones(Bubble.N,1);
%         for iBubble = 2:Bubble.N
%             q(1:iBubble-1) = sqrt(sum((Bubble.LocGlob(iBubble,:)- Bubble.LocGlob(1:iBubble-1,:)).^2,2));
%             
%             %         if (~sum(q==0)) ;break; end
%             Bubble.meandiff(iBubble) =  sum(q(1:iBubble-1)) ;
%         end
%         Bubble.avediff = sum(Bubble.meandiff)/(Bubble.N-1)*Bubble.N;
%     end
    %
    %     f = figure('WindowState','maximized');
    %     ax = gca; hold on;
    %     plot(domain.tpar ,Bubble.Contrast);
    %
    %     xlim([6e-6 12e-6])
    %     ax.YLabel.String = 'Pressure [Pa]';
    %     set(ax,'XTickLabel',ax.XTick*1e+6) ; ax.XLabel.String = 'Time [microsec]';
    %     title(['Comparison of Scattered Pressure of each microbubble' ])
    %
    %     disp(['All Bubbles have an inbetween distance higher than ', num2str(MinVal), ' [m]'])
    %
    %     %% Comparison of pressure initally and after converging
    %     disp('Calculating excitation pressure ...')
    %
    %     filename_pres = split(ls([file.dirname '/Bubbles/' 'output_file_pressure_*.txt']));
    %     filename_pres = {filename_pres{1:end-1}};
    %     filename_pres_init = filename_pres(length(filename_pres)/2+1:end);
    %     filename_pres_final = filename_pres(1:length(filename_pres)/2);
    %
    %     Bubble.InitExcPres = zeros(Bubble.N,domain.tdimpar);    Bubble. InitTime = zeros(Bubble.N,domain.tdimpar);
    %     Bubble.FinalExcPres = zeros(Bubble.N,domain.tdimpar);     Bubble.FinalTime = zeros(Bubble.N,domain.tdimpar);
    %
    %     for i  = 1:length(filename_pres)/2
    %         fileID_pres = fopen(filename_pres_init{i},'r');
    %         size_pres = [1 Inf];
    %         formatSpec = '%f';
    %         InitTimePres = fscanf(fileID_pres,formatSpec,size_pres);
    %         InitTimePres = reshape(InitTimePres,2*domain.tdimpar+1,length(InitTimePres)/(2*domain.tdimpar+1));
    %
    %         fileID_pres = fopen(filename_pres_final{i},'r');
    %         FinalTimePres = fscanf(fileID_pres,formatSpec,size_pres);
    %         FinalTimePres = reshape(FinalTimePres,2*domain.tdimpar+1,length(FinalTimePres)/(2*domain.tdimpar+1));
    %
    %         for BubbleNum = 1: size(InitTimePres,2)
    %             Bubble.InitExcPres(InitTimePres(1,BubbleNum),:) = InitTimePres(2:domain.tdimpar+1,BubbleNum);
    %             Bubble.InitTime(InitTimePres(1,BubbleNum),:) = InitTimePres(domain.tdimpar+2:end,BubbleNum);
    %         end
    %
    %         for BubbleNum = 1: size(FinalTimePres,2)
    %             Bubble.FinalExcPres(FinalTimePres(1,BubbleNum),:) = FinalTimePres(2:domain.tdimpar+1,BubbleNum);
    %             Bubble.FinalTime(FinalTimePres(1,BubbleNum),:) = FinalTimePres(domain.tdimpar+2:end,BubbleNum);
    %         end
    %     end
    %
end
disp(' ')
fclose('all')
end

function output = LoadContrastParams(Bubble,file,filename,input, load_operator, fformat)
if (nargin<5) 
    load_operator='no';
    fformat = 'r';
elseif (nargin<6)
    fformat = 'r';
end

if (strcmp(load_operator,'yes'))
    for iProcID =0:0
        
        filename_loc = [file.dirname '/Bubbles/',filename,int2string_ICS(iProcID)];
        fileID_loc = fopen(filename_loc,fformat);
        formatSpec = '%f %f %f';
        size_loc = [1 Inf];
        if strcmp(fformat,'r')  ;input(1+Bubble.N:2*Bubble.N,:) = reshape(fscanf(fileID_loc,formatSpec,size_loc),size(input,2),Bubble.N)'; end
        if strcmp(fformat,'rb') ;input(1+Bubble.N:2*Bubble.N,:) = reshape(fread(fileID_loc,size_loc,'*double')  ,size(input,2),Bubble.N)'; end
        %         Check if all the values of all the processors are the same
        if iProcID >0
            if sum(sum(abs(input(1+Bubble.N:2*Bubble.N,:)-input(1:Bubble.N,:))>1E-8))~=0 ;error(['ERROR in ',filename]); return; end
        end
        input(1:Bubble.N,:) = input(1+Bubble.N:2*Bubble.N,:) ;
    end
end
output = input(1:Bubble.N,:);
end

function output = LoadContrastParams_F90_BIN(Bubble,file,filename,input, load_operator)
% File arrtimes.m
% x = fread(fid,[10 1],'*double')

if (strcmp(load_operator,'yes'))
    for iProcID =0:0
        
        filename_loc = [file.dirname '/Bubbles/',filename,int2string_ICS(iProcID)];
        fileID_loc = fopen(filename_loc,'rb');
%         DISCARD_FIRST_BYTES = fread(fileID_loc,[1 1],'*int64');
        formatSpec = '%f %f %f';
        size_loc = [1 Inf];
        input(1+Bubble.N:2*Bubble.N,:) = reshape(fread(fileID_loc,size_loc,'*double'),size(input,2),Bubble.N)';
        %         Check if all the values of all the processors are the same
        if iProcID >0
            if sum(sum(abs(input(1+Bubble.N:2*Bubble.N,:)-input(1:Bubble.N,:))>1E-8))~=0 ;error(['ERROR in ',filename]); return; end
        end
        input(1:Bubble.N,:) = input(1+Bubble.N:2*Bubble.N,:) ;
    end
end
output = input(1:Bubble.N,:);
end