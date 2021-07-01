function [lt_opt,thetat_opt,st_opt]=check_ICS_timeframe(mediumparams,domainparams,...
    modelparams,sourceparams,showplot)

% [lt_opt,thetat_opt,st_opt]=check_ICS_timeframe(mediumparams,domainparams,
%                                                  modelparams,sourceparams,showplot)
% Check whether the defined time frame will fit the field excitation. The
% first and last arrival times at the upper and lower domain boundaries 
% in the x-direction and the y-direction of the computational domain are
% checked. There are some small simplifications that may result in slightly 
% different actual first and last arrival times, but these will be small.
% Inputs:
% mediumparams: structure that has been generated with set_ICS_medium.
% domainparams: structure that has been generated with set_ICS_domain.
% modelparams:  structure that has been generated with set_ICS_model. 
% sourceparams: structure that has been generated with set_ICS_source_parametric.
%               Check_ICS_timeframe only works for parametric aperture and 
%               signature definitions, not for data/file definitions.
% showplot:     (optional, default = 1) whether (1) or not (0) to plot the
%               current and optimal timeframe.
% Outputs:
% lt_opt:       optimal time length that makes the full excitation
%               fit in the co-moving time frame. [1/freq0]
% thetat_opt:      optimal angle of the time frame [rad]
% st_opt:       optimal start of the excitation [1/freq0]

% JH - Version 090121

% check inputs
if(strcmp(sourceparams.definitiontype,'file'))
    disp('Warning: check_ICS_timeframe only works for parametric source definitions')
    return
elseif(strcmp(sourceparams.signature.signaturetype,'file'))
    disp('Warning: check_ICS_timeframe only works for parametric source signature definitions')
    return
end

% Be sure to get all time parameters normalized with respect to 1/freq0
[lt,st] = normalize_lt_st_to_freq0(domainparams,modelparams.freq0);

% Normalize f0 to be counted in MHz, Tpulse to be counted in microseconds, 
% and c0 to be counted in mm/mus
freq0 =  modelparams.freq0/1e6;
c0 =     mediumparams.c0*1e3/1e6;
Tpulse = sourceparams.signature.Tpulse/freq0;
lt =     lt/freq0;
st =     st/freq0;

% correct lt to include st, as in generate_ICS_configurationfile()
lt = lt + st;

% Obtain the z-axis and obtain the x-axis as a function of the z-axis
dt =    1/(2*modelparams.Fnyq*freq0);
dx =    c0*dt;      %all distance parameters are in mm
zdim =  ceil(domainparams.lz/dx);
zaxis = (1:zdim)*dx;
xaxis = zaxis/tan(domainparams.thetax);
yaxis = zaxis/tan(domainparams.thetay);

% shorten notation
ap = sourceparams.aperture;
lx = domainparams.lx;
ly = domainparams.ly;

if (strcmp(sourceparams.aperture.aperturetype,'cylindrical'))
    
    % Calculate first and last time of arrival from five points on 
    % the cylindrical aperture (center and extremes in x and y)
    % to each z-position on the four domain boundaries, then determine min and max.
    % Simplification: only for five points, not for all points on the
    % aperture
    
    xsrc = [-1 1  0 0 0] * ap.radius;
    ysrc = [ 0 0 -1 1 0] * ap.radius;
    %account for positive, negative or no focus
    if(abs(ap.radfocus)>1e-5)
        if(ap.radfocus>0)
            tdelaysrc = [0 0 0 0 1]*(sqrt(ap.radius^2+ap.radfocus^2)-abs(ap.radfocus))/c0;
        else
            tdelaysrc = [1 1 1 1 0]*(sqrt(ap.radius^2+ap.radfocus^2)-abs(ap.radfocus))/c0;
        end
    else
        tdelaysrc = [0 0 0 0 0];
    end

    [tdelaysrcmat,xaxismat] = meshgrid(tdelaysrc,xaxis);
    [tdelaysrcmat,yaxismat] = meshgrid(tdelaysrc,yaxis);
    [xsrcmat,zaxismat] =      meshgrid(xsrc,zaxis);
    [ysrcmat,zaxismat] =      meshgrid(ysrc,zaxis);

    % firstTOA takes first time of arrival from all source points to all z 
    % positions over the beam on the upper and lower boundaries in the x 
    % direction (first and second), and in the upper and lower boundaries 
    % in the y direction (third and fourth)
    firstTOA = min([...
        min(tdelaysrcmat+sqrt((xaxismat-lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(tdelaysrcmat+sqrt((xaxismat+lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(tdelaysrcmat+sqrt((xaxismat-xsrcmat).^2+(yaxismat-ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(tdelaysrcmat+sqrt((xaxismat-xsrcmat).^2+(yaxismat+ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'],[],1);
    % lastTOA takes last time of arrival etc. etc..
    lastTOA = Tpulse + max([...
        max(tdelaysrcmat+sqrt((xaxismat-lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(tdelaysrcmat+sqrt((xaxismat+lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(tdelaysrcmat+sqrt((xaxismat-xsrcmat).^2+(yaxismat-ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(tdelaysrcmat+sqrt((xaxismat-xsrcmat).^2+(yaxismat+ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'],[],1);

elseif (strcmp(sourceparams.aperture.aperturetype,'rectangular'))
    
    % Calculate first and last time of arrival from four points on 
    % the rectangular aperture (extremes in x and y)
    % to each z-position on the four domain boundaries, then determine min and max.
    % Simplification: only for four points, not for all points on the
    % aperture
    
    xsrc = [-1 1  0 0] * ap.rectwidth/2;
    ysrc = [ 0 0 -1 1] * ap.rectheight/2;

    [xsrcmat,xaxismat] = meshgrid(xsrc,xaxis);
    [ysrcmat,yaxismat] = meshgrid(ysrc,yaxis);
    [xsrcmat,zaxismat] = meshgrid(xsrc,zaxis);

    % firstTOA takes first time of arrival from all source points to all z 
    % positions over the beam on the upper and lower boundaries in the x 
    % direction (first and second), and in the upper and lower boundaries 
    % in the y direction (third and fourth)
    firstTOA = min([...
        min(sqrt((xaxismat-lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(sqrt((xaxismat+lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(sqrt((xaxismat-xsrcmat).^2+(yaxismat-ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(sqrt((xaxismat-xsrcmat).^2+(yaxismat+ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'],[],1);
    % lastTOA takes last time of arrival etc. etc..
    lastTOA = Tpulse + max([...
        max(sqrt((xaxismat-lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(sqrt((xaxismat+lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(sqrt((xaxismat-xsrcmat).^2+(yaxismat-ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(sqrt((xaxismat-xsrcmat).^2+(yaxismat+ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'],[],1);

elseif (strcmp(sourceparams.aperture.aperturetype,'triangular'))
    
    % Calculate first and last time of arrival from four points on 
    % a rectangle around the triangular aperture (extremes in x and y)
    % to each z-position on the four domain boundaries, then determine min and max.
    % Simplification: treat it as though it is a rectangular-shaped source,
    % check only for four points, not for all points on the aperture
    
    xsrc = [-1 1  0 0] * sqrt(3)/2*ap.triangedge + ap.triangxcenter;
    ysrc = [ 0 0 -1 1] * ap.triangedge/2;

    [xsrcmat,xaxismat] = meshgrid(xsrc,xaxis);
    [ysrcmat,yaxismat] = meshgrid(ysrc,yaxis);
    [xsrcmat,zaxismat] = meshgrid(xsrc,zaxis);

    % firstTOA takes first time of arrival from all source points to all z 
    % positions over the beam on the upper and lower boundaries in the x 
    % direction (first and second), and in the upper and lower boundaries 
    % in the y direction (third and fourth)
    firstTOA = min([...
        min(sqrt((xaxismat-lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(sqrt((xaxismat+lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(sqrt((xaxismat-xsrcmat).^2+(yaxismat-ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        min(sqrt((xaxismat-xsrcmat).^2+(yaxismat+ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'],[],1);
    % lastTOA takes last time of arrival etc. etc..
    lastTOA = Tpulse + max([...
        max(sqrt((xaxismat-lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(sqrt((xaxismat+lx/2-xsrcmat).^2+(yaxismat-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(sqrt((xaxismat-xsrcmat).^2+(yaxismat-ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'; ...
        max(sqrt((xaxismat-xsrcmat).^2+(yaxismat+ly/2-ysrcmat).^2+zaxismat.^2)/c0,[],2)'],[],1);

elseif (strcmp(sourceparams.aperture.aperturetype,'pointsource'))
    
    % Calculate first and last time of arrival from the point source
    % in (x,y) = (0,0) to each z-position on the four domain boundaries, 
    % then determine min and max.
    
    % firstTOA takes first time of arrival from all source points to all z 
    % positions over the beam on the upper and lower boundaries in the x 
    % direction (first and second), and in the upper and lower boundaries 
    % in the y direction (third and fourth)
    firstTOA = min([...
        sqrt((xaxis-lx/2).^2+(yaxis).^2+zaxis.^2)/c0; ...
        sqrt((xaxis+lx/2).^2+(yaxis).^2+zaxis.^2)/c0; ...
        sqrt((xaxis).^2+(yaxis-ly/2).^2+zaxis.^2)/c0; ...
        sqrt((xaxis).^2+(yaxis+ly/2).^2+zaxis.^2)/c0],[],1);
    % lastTOA takes last time of arrival etc. etc..
    lastTOA = Tpulse + firstTOA;
    
elseif (strcmp(ap.aperturetype,'phasedarray'))

    % Calculate first and last time of arrival for each array element to each
    % z-position on the domain boundaries, then determine min and max.
    % Simplifications: on the upper and lower boundaries in the x-direction 
    % only from the center of each element, and neglecting the elevation focus

    % obtain the positions and phase delays of the source elements
    half_arraywidth = .5*(ap.numels*(ap.elwidth+ap.kerf)-ap.kerf);
    elx =             -half_arraywidth+ap.elwidth/2+(0:ap.numels-1)*(ap.elwidth+ap.kerf);
    eltdelay =        -sqrt((ap.focusx-elx).^2+ap.focusz^2);
    eltdelay =        eltdelay-min(eltdelay);

    [eltdelaymat,xaxismat] = meshgrid(eltdelay,xaxis);
    [eltdelaymat,yaxismat] = meshgrid(eltdelay,yaxis);
    [elxmat,zaxismat] =      meshgrid(elx,zaxis);

    % firstTOA takes first time of arrival from all elements to all z 
    % positions over the beam on the upper and lower boundaries in the x 
    % direction (first and second), and in the upper and lower boundaries 
    % in the y direction (third and fourth)
    firstTOA = min([...
        min(eltdelaymat+sqrt((xaxismat-lx/2-elxmat).^2+yaxismat.^2+zaxismat.^2)/c0,[],2)'; ...
        min(eltdelaymat+sqrt((xaxismat+lx/2-elxmat).^2+yaxismat.^2+zaxismat.^2)/c0,[],2)'; ...
        min(eltdelaymat+sqrt((xaxismat-elxmat).^2+(yaxismat-ly/2-ap.elheight/2).^2+zaxismat.^2)/c0,[],2)'; ...
        min(eltdelaymat+sqrt((xaxismat-elxmat).^2+(yaxismat+ly/2+ap.elheight/2).^2+zaxismat.^2)/c0,[],2)'],[],1);
    % lastTOA takes last time of arrival etc. etc..
    lastTOA = Tpulse + max([...
        max(eltdelaymat+sqrt((xaxismat-lx/2-elxmat).^2+yaxismat.^2+zaxismat.^2)/c0,[],2)'; ...
        max(eltdelaymat+sqrt((xaxismat+lx/2-elxmat).^2+yaxismat.^2+zaxismat.^2)/c0,[],2)'; ...
        max(eltdelaymat+sqrt((xaxismat-elxmat).^2+(yaxismat-ly/2-ap.elheight/2).^2+zaxismat.^2)/c0,[],2)'; ...
        max(eltdelaymat+sqrt((xaxismat-elxmat).^2+(yaxismat+ly/2+ap.elheight/2).^2+zaxismat.^2)/c0,[],2)'],[],1);

end

% Calculate the T beam as in the domainparams structure
deltat = zaxis/c0/tan(domainparams.thetat) - st;

% Calculate optimized T beam - calculate angle as though the pulse comes 
% from the center of the transducer - is correct for the far field
thetat_opt = atan(1/sqrt(1/cos(pi/2-domainparams.thetax)^2+tan(pi/2-domainparams.thetay)^2));
deltat_opt = zaxis/c0/tan(thetat_opt);
st_opt= - min([ (firstTOA-deltat_opt) 0]);
deltat_opt = deltat_opt - st_opt;
lt_opt = max(lastTOA-deltat_opt);

%Plot figure of the beams and arrival times
if(showplot~=0)
    figure(10),plot(...
        zaxis,firstTOA*freq0,'b',...
        zaxis,lastTOA*freq0,'b-.',...
        zaxis,deltat*freq0,'g',...
        zaxis,deltat_opt*freq0,'c:',...
        zaxis,(deltat+lt)*freq0,'g-.',...
        zaxis,(deltat_opt+lt_opt)*freq0,'c:')
    xlabel('z (mm)')
    ylabel('t (periods re f_0)')
    legend('First TOA','Last TOA','Current T beam','Optimal T beam','Location','SouthEast')
end

% set output variables to correct normalization
lt_opt =  lt_opt*freq0;
thetat_opt = thetat_opt;
st_opt =  st_opt*freq0;
