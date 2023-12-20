function [maxhilb,output]=maxhilbert_vec(input,dt,lowfreq,highfreq,order,isperiodical, interpval)

if(nargin<4)
    error('Not enough arguments')
elseif(nargin<5)
    order=8;
end
if(nargin<6)
    isperiodical=false;
end
if (nargin<7)
    interpval = 4;
elseif(nargin>7)
    error('Too much arguments')
end

dims=size(input);

maxfreq=1/(2*dt);

Wn = [lowfreq highfreq]/maxfreq;
[B,A]=butter(order,[lowfreq highfreq]/maxfreq);
[A1,B1,C1,D1]=butter(order,Wn);

if(isperiodical)
    initiallength=max(length(A),length(B))-1;
    initialcond=input((-initiallength+1:0)+dims(1),:);
    output=filter(B,A,input,initialcond);
    
    output=statefilter(input,A1,B1,C1,D1);
    output=output(initiallength+(1:dims(1)),:);
else
    output=filter(B,A,input);  % filters the input data using a rational transfer function defined by the numerator and denominator coefficients b and a.
end

if (interpval>1)
    output=interpft(output,interpval*size(output,1));
end
maxhilb=squeeze(max(abs(hilbert(output))));

% if (length(dims)>2)
%     output=reshape(output,dims);
%     maxhilb=reshape(maxhilb,dims(2:length(dims)));
% end
% keyboard

