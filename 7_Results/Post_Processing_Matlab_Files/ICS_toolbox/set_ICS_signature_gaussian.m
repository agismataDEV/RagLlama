function signature = set_ICS_signature_gaussian(Pstart,Tpulse,Twidth,Tdelay,...
    power,freq0)

% signature = set_ICS_signature_gaussian(Pstart,Tpulse,Twidth,Tdelay,power,freq0)
% Generates the signature structure for a harmonic signal modulated 
% with a gaussian window. The signature structure contains the 
% signature settings used in the set_ICS_source_parametric function. 
% Inputs:
% Pstart: source pressure amplitude in [Pa].
% Tpulse: total length of the pulse in [1/freq0].
% Twidth: width of the Gaussian window in [1/freq0].
% Tdelay: delay of the Gaussian window in [1/freq0].
% power:  (optional, default = 2) power in the exponent of the Gaussian window.
% freq0:  (optional, default = 0) center frequency of the harmonic 
%         signal in [Hz]. If freq0=0, then the center frequency of the 
%         harmonic signal is then set to the fundamental frequency set in 
%         set_ICS_model.

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(nargin<5 || isempty(power))
    power = 2;
end
if(nargin<6 || isempty(freq0))
    freq0 = 0;
end
if(nargin>6)
    error('Too many arguments')
end

if((Pstart<=0)||(Tpulse<=0)||(Twidth<=0)||(power<=0)||(freq0<0))
    error('Invalid value for any of the parameters')
end

signature = struct('signaturetype','gaussian','freq0',freq0,...
    'Pstart',Pstart,'Tpulse',Tpulse,'Twidth',Twidth,...
    'Tdelay',Tdelay,'power',power);