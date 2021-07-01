function [lt,st] = normalize_lt_st_to_freq0(domainparams,freq0)

% [lt,st] = normalize_lt_st_to_freq0(domainparams,modelparams)
% check freq0 in domainparams, and if
% necessary, normalize lt and st to 1/freq0 instead of microseconds

% JH - version 090121

if(domainparams.freq0~=0)
    if(domainparams.freq0~=freq0)
        error('freq0 in signature of domainparams is incompatible with freq0 of modelparams')
    end
    lt = domainparams.lt;
    st = domainparams.st;
else
    lt = domainparams.lt*1e-6*freq0;
    st = domainparams.st*1e-6*freq0;
end
