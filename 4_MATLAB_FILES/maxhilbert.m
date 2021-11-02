function [maxhilb,output]=maxhilbert(input,dt,lowfreq,highfreq,order,isperiodical)

%[maxhilb,output]=maxhilbert(input,dt,lowfreq,highfreq,order,isperiodical)

if(nargin<4)
    error('Not enough arguments')
elseif(nargin<5)
    order=8;
end
if(nargin<6)
    isperiodical=false;
elseif(nargin>6)
    error('Too much arguments')
end   

dims=size(input);
input=input(:,:);

maxfreq=1/(2*dt);

[B,A]=butter(order,[lowfreq highfreq]/maxfreq);

if(isperiodical)
    initiallength=max(length(A),length(B))-1;
    initialcond=input((-initiallength+1:0)+dims(1),:);
    output=filter(B,A,input,initialcond);
else
    output=filter(B,A,input);
end

maxhilb=squeeze(max(abs(hilbert(output))));

if (length(dims)>2)
    output=reshape(output,dims);
    maxhilb=reshape(maxhilb,dims(2:length(dims)));
end
% keyboard

