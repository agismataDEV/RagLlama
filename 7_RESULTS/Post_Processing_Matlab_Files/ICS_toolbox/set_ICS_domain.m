function domainparams = set_ICS_domain(lt,lx,ly,lz,st,sz,thetat,thetaxy,freq0)

% domainparams = set_ICS_domain(lt,lx,ly,lz,st,sz,thetat,thetaxy,freq0)
% Generates the domainparams structure that contains all information on the
% domain size and shape used in the ICS Parnac program. 
% Inputs:
% lt:      length of the time axis [1/freq0] or [\mus]
% lx:      length of the x axis [mm]
% ly:      length of the y axis [mm]
% lz:      length of the z axis [mm]
% st:      (optional, default = 0) start time of the excitation [1/freq0] or [\mus]
%          beware: the time frame itself always starts at t=0. If st>0, lt
%          will be corrected in generate_ICS_configurationfile to lt = lt+st.
% sz:      (optional, default = 0) start position in z of the domain [mm]
% thetat:  (optional, default = pi/4) angle in radians, in the range
%          [0 pi/2], of the computational domain in the t-z plane,
%          measured with respect to the positive t axis [rad]
% thetaxy: (optional, default = [pi/2 pi/2]) vector of two elements
%          [thetax thetay] denoting the angles in radians, both in the range
%          [0 pi], of the projection of the computational domain on the 
%          x-z and y-z plane, both measured with respect to the positive 
%          x and y axes
% freq0:   (optional, default = 0) fundamental frequency to which 
%          lt and st may be related. If freq0=0, then lt and st are
%          given in terms of microseconds instead of in 1/freq0.

% JH, version 090116

if(nargin<4)
    error('Not enough arguments')
end
if(nargin<5 || isempty(st))
    st = 0;
end
if(nargin<6 || isempty(sz))
    sz = 0;
end
if(nargin<7 || isempty(thetat))
    thetat = pi/4;
end
if(nargin<8 || isempty(thetaxy))
    thetaxy = [pi/2 pi/2];
end
if numel(thetaxy)~=2
    error('Thetaxy should be a vector of two elements')
end
if(nargin<9|| isempty(freq0))
    freq0 = 0;
end
if(nargin>9)
    error('Too many arguments')
end

if((lt<=0)||(lx<=0)||(ly<=0)||(lz<=0)||(st<0)||(sz<0))
    error('One of the parameters has a negative or otherwise invalid value')
end
if(thetat<0 || thetat>pi/2)
    error('Thetat should be between 0 and pi/2')
end
if(min(thetaxy)<0 || max(thetaxy)>pi)
    error('Both elements of thetaxy should be between 0 and pi')
end

domainparams = struct('lt',lt,'lx',lx, 'ly', ly, 'lz', lz, ...
    'st', st, 'sz', sz, 'thetat', thetat,...
    'thetax', thetaxy(1), 'thetay', thetaxy(2),'freq0',freq0);