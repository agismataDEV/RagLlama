function memsize = get_ICS_arraysizes(mediumparams,domainparams,...
    modelparams,sourceparams,numprocs,dTdelay,iDimSource)

% memsize = get_ICS_arraysizes(mediumparams,domainparams,...
%     modelparams,sourceparams,numprocs,dTdelay,iDimSource)
% Calculates the arraysizes used internally in the ICS Parnac program

st,sx,sy,lt,lx,ly,lz,tht,thx,thy,...
    sourcetype,FDXorder,f0,fnyq,c0,numproc)

% memsize =
% critical_memorysize_withsource(st,sx,sy,lt,lx,ly,lz,tht,thx,thy,sourcetype,FDXorder,f0,fnyq,c0,numproc)
% calculates the critical memorysize for the ICS Parnac program

if(abs(thy-pi/2)<1e-5)
    ly=ly/2;
end

%normalizes and translates to number of points
[ist,isx,isy,it,ix,iy,iz]=translate_to_points_withsource(st,sx,sy,lt,lx,ly,lz,f0,fnyq,c0);

%size of a field, resized to 'nice' values in terms of prime numbers (FFT) 
% and number of processors
disp('Near-field and beam space old and new sizes:')
disp([' ' num2str(it) ' ' num2str(ix) ' ' num2str(iy) ' ' num2str(iz)]);
[spt,spx,spy,spz]=correct_size(it,ix,iy,iz,numproc);
disp([' ' num2str(spt) ' ' num2str(spx) ' ' num2str(spy) ' ' num2str(spz)]);

%size of a field including dp/dz margins
tant=1/tan(tht);
tanx=1/tan(thx);
tany=1/tan(thy);

%extra margins for linear field solution is only needed in case of an
%F-source
if(sourcetype~=1)
    FDXorder = 0;
end

ist = max(ist,it)+ 2*ceil(abs(FDXorder/2*tant)-1e-5);

disp('Source space old and new sizes:')
disp([' ' num2str(ist) ' ' num2str(isx) ' ' num2str(isy) ' 1' ]);
[sspt,sspx,sspy,sspz]=correct_size(ist,isx,isy,1,numproc);
disp([' ' num2str(sspt) ' ' num2str(sspx) ' ' num2str(sspy) ' ' num2str(sspz)]);

mt = ist;
mx = spx+2*ceil(abs(FDXorder/2*tanx)-1e-5);
my = spy+2*ceil(abs(FDXorder/2*tany)-1e-5);
mz = spz+2*FDXorder/2;

disp('Nearfield/beam space with TXYZ margin old and new sizes:')
disp([' ' num2str(mt) ' ' num2str(mx) ' ' num2str(my) ' ' num2str(mz)]);
[mt,mx,my,mz]=correct_size(mt,mx,my,mz,numproc);
disp([' ' num2str(mt) ' ' num2str(mx) ' ' num2str(my) ' ' num2str(mz)]);

%size of slices, with and without margin
%assumptions:
% - source and dest slice have the same angle, therefore
%       no extra skew margins are accounted for
% - source and dest slice have the same size for 
%       field->field calculations
% primary source: is accounted for correctly in the _withsource form

%without margins
slx = 2 * spx - 1;
sly = 2 * spy - 1;
slz = 2 * spz - 1; 
while (~Isprodsmallprimes(slx))
        slx=slx+1;
end
while (~Isprodsmallprimes(sly))
        sly=sly+1;
end
while (~Isprodsmallprimes(slz))
        slz=slz+1;
end

disp('One nearfield/beam slice old and new sizes:')
disp([' 1 ' num2str(slx) ' ' num2str(sly) ' ' num2str(slz)]);
%[temp,slx,sly,slz]=correct_size(spt,slx,sly,slz,numproc);
disp([' 1 ' num2str(slx) ' ' num2str(sly) ' ' num2str(slz)]);

%with margins
slmx = mx + sspx - 1;
slmy = my + sspy - 1;
slmz = mz; %since this is only used to calculate the linear field solution

while (~Isprodsmallprimes(slmx))
        slmx=slmx+1;
end
while (~Isprodsmallprimes(slmy))
        slmy=slmy+1;
end
while (~Isprodsmallprimes(slmz))
        slmz=slmz+1;
end

disp('One Nearfield/beam slice with TXYZ margin old and new sizes:')
disp([' 1 ' num2str(slmx) ' ' num2str(slmy) ' ' num2str(slmz)]);
%[temp,slmx,slmy,slmz]=correct_size(spt,slmx,slmy,slmz,numproc);
disp([' 1 ' num2str(slmx) ' ' num2str(slmy) ' ' num2str(slmz)]);

%total memory taken on each processor by the different arrays, in Megabytes
sourcedist12size = sspt*sspx*sspy*sspz*2*8/1024^2/numproc
fielddist0size = spt*spx*spy*spz*8/1024^2/numproc
fielddist12size = (spt+1)*spx*spy*spz*2*8/1024^2/numproc
fieldwithmargindist12size = (mt+1)*mx*my*mz*2*8/1024^2/numproc
slicewithmargindist2size = slmx*slmy*slmz*2*8/1024^2
slicedist2size = slx*sly*slz*2*8/1024^2

%maximum used memory in Megabytes on specific instances in the program
% 1st: in the linear NF solution, during the convolution (neglect smaller
% slice...)
% 2nd: in the linear NF solution, during the 2->1 redist. afterwards
% 3rd: in the nonlinear NF solution, during the convolution

memsize_nobeam = [ ...
    sourcedist12size + fieldwithmargindist12size + slicewithmargindist2size ...
    sourcedist12size + 2*fieldwithmargindist12size ...
    fielddist12size + 2*slicedist2size ...
    ]

%maximum used memory in Megabytes on specific instances in the program
% 1st: in the nonlinear NF or BEAM solution, during the convolution
% 2nd: in the linear BEAM solution, during the convolution
% 3rd: in the linear BEAM solution, during the 2->1 redist. afterwards
% 4th: in the NF->BEAM interaction, during the convolution
% 5th: in the NF->BEAM interaction, during the 1->2 redist.

memsize = [ ...
    fielddist12size + 2*slicedist2size ...
    fieldwithmargindist12size + slicewithmargindist2size + fielddist0size ...
    2*fieldwithmargindist12size + fielddist0size ...
    2*fielddist12size + 2*slicedist2size ...
    3*fielddist12size ...
    ]

%%get_ICS_discretesizes

function [ist,isx,isy,it,ix,iy,iz]=get_ICS_discretesizes(st,sx,sy,lt,lx,ly,lz,f0,fnyq,c0)

% [ist,isx,isy,it,ix,iy,iz]=get_ICS_discretesizes(st,sx,sy,lt,lx,ly,lz,f0,fnyq,c0)

%to normalized form
snt = st;
snx = sx*1e-3*f0/c0;
sny = sy*1e-3*f0/c0;
lnt = lt;
lnx = lx*1e-3*f0/c0;
lny = ly*1e-3*f0/c0;
lnz = lz*1e-3*f0/c0;

%to number of points
dt = 1/(2*fnyq);
dx = 1/(2*fnyq);

ist = ceil(snt/dt);
isx = floor(2*(snx/(2*dx) + 10 ) + 1); 
isy = floor(2*(sny/(2*dx) + 10 ) + 1); 
% for phased array, we would need + 15 instead of + 10...
it = ceil(lnt/dt - 1e-5);
ix = ceil(lnx/dx - 1e-5); 
iy = ceil(lny/dx - 1e-5); 
iz = ceil(lnz/dx - 1e-5) - 1;

%% end get_ICS_discretesizes

%% get_ICS_correctedsizes

function [pt,px,py,pz]=get_ICS_correctedsizes(it,ix,iy,iz,numproc)

% [pt,px,py,pz]=correct_size(it,ix,iy,iz,numproc)
% implements the f90 function InitGridSpace

%resize to 'nice' number of points (for numprocs)
itnew = ceil((it+1)/numproc)*numproc - 1;

if (Isprodsmallprimes(ix*iy*iz) & (mod(ix*iy*iz,numproc)==0))
    ixnew = ix;
    iynew = iy;
    iznew = iz;
else
    
    iPa                = 1;
    iPb                = 1;
    iPc                = 1;
    i                = numproc;
    j                = 2; %first prime number
    k                = 0; %counter for the xyz-dimensions
    while (i ~= 1)
        % If not j | i, we should try the next number
        while (mod(i, j) ~= 0)
            j = j + 1;
        end

        if (mod(k,3)==0)
            iPa = iPa * j;
        elseif (mod(k,3)==1)
            iPb = iPb * j;
        elseif (mod(k,3)==2)
            iPc = iPc * j;
        end

        i        = floor(i / j);
        k         = k + 1;
    end
    
    divx = ceil(ix/iPa);
    while (~Isprodsmallprimes(divx))
        divx=divx+1;
    end
    divy = ceil(iy/iPb);
    while (~Isprodsmallprimes(divy))
        divy=divy+1;
    end
    divz = ceil(iz/iPc);
    while (~Isprodsmallprimes(divz))
        divz=divz+1;
    end

    ixnew = divx * iPa;
    iynew = divy * iPb;
    iznew = divz * iPc;

end

pt = itnew;
px = ixnew;
py = iynew;
pz = iznew;

%%end get_ICS_correctedsizes
                
        