% This function is used mainly for fixing the spectral leakage in the lower harmonics
% It's a hard cutoff filter in the frequency domain. A blackman window is applied
% in the time domain. With this method, the scattered pressure of one small scatterer
% can be visualized in the lower harmonics.

function output=Bubble_FilterINCS(input,dt,PPW_t,lowfreq,highfreq,order,isperiodical)

if(nargin<5)
    error('Not enough arguments')
elseif(nargin<6)
    order=8;
end
if(nargin<7)
    isperiodical=false;
elseif(nargin>7)
    error('Too much arguments')
end   

dims=size(input);
w  = ones(dims(1),dims(2));
w  = repmat(blackman(dims(1)),1,dims(2));
input=input(:,:).*w;

maxfreq = 1/(PPW_t*dt);
n=2^nextpow2(size(input,1));
f = maxfreq*(0:n/2)/(n/2); 
input = abs(fft(input,n));
input =input(1:n/2+1,:);

lowbound = lowfreq;
highbound = highfreq;
% % w  = ones(size(find(f>lowfreq & f<=highfreq),2),size(input,2));
% w  = repmat(hann(size(input,1)),1,size(input,2)); % Hann window in frequency
% w  = repmat(hann(size(input,2)),1,size(input,1))'; % Hann window in space
w  = repmat(hann(size(find(f>lowbound & f<=highbound),2)),1,size(input,2)); % Hann window in frequency
input(find(f>=lowbound & f<=highbound),:) = input(find(f>=lowbound & f<=highbound),:).*w;
input(find(f<lowbound | f>highbound),:) = 0;
input=interpft(input,1*size(input,1));

output = squeeze(max(abs(ifft(input))));
end 
