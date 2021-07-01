function signature = set_ICS_signature_blackman(Pstart,Tpulse, freq0)

% signature = set_ICS_signature_blackman(Pstart,Tpulse)
% Generates the signature structure for a blackman pulse. 
% The signature structure contains the signature settings 
% used in the set_ICS_source_parametric function. 
% Inputs:
% Pstart: source pressure amplitude in [Pa].
% Tpulse: total length of the pulse in [1/freq0] or [\mus].
% freq0:  (optional, default = 0) center frequency of the harmonic 
%         signal in [Hz]. If freq0=0, then Tpulse is
%         given in terms of microseconds instead of in 1/freq0.

% JH, version 090116

if(nargin<1)
    error('Not enough arguments')
end
if(nargin<3 || isempty(freq0))
    freq0 = 0;
end
if(nargin>3)
    error('Too many arguments')
end

if((Pstart<=0)||(Tpulse<=0)||(freq0<=0))
    error('Invalid value for any of the parameters')
end

signature = struct('signaturetype','blackman','freq0',freq0,...
    'Pstart',Pstart,'Tpulse',Tpulse,'Tdelay',0);