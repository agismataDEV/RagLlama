%% L.D. 10-05-2012 - This script reads and plots the residuals
clear all
close all
clc
concentrated = '_concentrated';
% concentrated='';
% plots
figure(1)
BubbleN = [1,10,100,1000];
% hold on
for Bubblei = 1:length(BubbleN)
name=[int2str(BubbleN(Bubblei)),'_qsource',concentrated,'\Neumann'] ;                  % file name
% name='BiCGSTAB'                  % file name
% name='SteepestDescent'           % file name

ProcessorsNumber=8;              % number of processors used for the simulation

for i=0:ProcessorsNumber-1
    if i<10
    err = textread(['.\',name,'errorN000000000',num2str(i),'.dat']); 
    else
    err = textread(['.\',name,'errorN00000000',num2str(i),'.dat']);  
    end
    Dim=size(err);
    TotalError(:,i+1)=err(:,1);
end

error=sum(TotalError');

for i=1:Dim(1)-1
    Residual(i+1)=sqrt(error(i+1)./error(1));
end
Residual(1)=1;

i_end = 's';
if Bubblei==1 ; i_end='';end
legend_data1{Bubblei} = [int2str(BubbleN(Bubblei)), ' Bubble',i_end ' /1st Iter'];
Error1(Bubblei,:) = Residual;
end


% Create multiple lines using matrix input to plot
plot1 = semilogy(0:Dim(1)-1,Error1);
% axes1 = axes;
grid on;
set(plot1(1),'DisplayName','1 Bubble','Marker','o','LineWidth',1,...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
set(plot1(2),'DisplayName','10 Bubbles','Marker','x',...
     'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
set(plot1(3),'DisplayName','100 Bubbles','Marker','square',...
     'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);
set(plot1(4),'DisplayName','1000 Bubbles','Marker','diamond',...
     'Color',[0.149019607843137 0.149019607843137 0.149019607843137]);

% Create ylabel
ylabel('residual log[Err(n)]');

% Create xlabel
xlabel('Number of Iterations [n]');
% Set the remaining axes properties
% legend(legend_data1)

%% ERR1(n) analisys
% Create axes

% figure(2)
hold on;


for Bubblei = 1:length(BubbleN)
name=[int2str(BubbleN(Bubblei)),'_qsource',concentrated,'\TESTNeumann'] ;  
first  = [name '_000_y_001_000'];
output = load_parnac_slice(first); 
plin   = squeeze(output.data); 
IterationNumber = 9;
ii = 0;
for jj=1:IterationNumber
    
    
    third       = [name '_00' num2str(jj) '_y_001_00' num2str(ii)];

    if (jj>9)&&(jj<100)
            third       = [name '_0' num2str(jj) '_y_001_00' num2str(ii)];
    else if jj>99
            third       = [name '_' num2str(jj) '_y_001_00' num2str(ii)];
        end
    end
    
    output = load_parnac_slice(third); 
    pnonl  = squeeze(output.data); 
    
    first       = [name '_00' num2str(jj-1) '_y_001_00' num2str(ii)];

    if (jj-1>9)&&(jj-1<100)
            first       = [name '_0' num2str(jj-1) '_y_001_00' num2str(ii)];
    else if jj-1>99
            first       = [name '_' num2str(jj-1) '_y_001_00' num2str(ii)];
        end
    end
        
    output = load_parnac_slice(first); 
    pnonl2 = squeeze(output.data); 
    
    Error2(Bubblei,jj) = sqrt(sum(sum(sum((pnonl2-pnonl).^2)))./sum(sum(sum((plin).^2))));
end

% PLOTS
i_end = 's';
if Bubblei==1 ; i_end='';end
legend_data2{Bubblei} = [int2str(BubbleN(Bubblei)), ' Bubble',i_end ' /Prev Iter'];
%  plot(0:IterationNumber,log([1,Error]))
Error2(Bubblei,1:10) = ([1,Error2(Bubblei,1:9)]);

end

grid on;

% Create multiple lines using matrix input to plot
plot2 = semilogy(0:IterationNumber,Error2);
set(plot2(1),'DisplayName','1 Bubble','Marker','o','LineWidth',1,...
    'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
set(plot2(2),'DisplayName','10 Bubbles','Marker','x',...
    'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
set(plot2(3),'DisplayName','100 Bubbles','Marker','square',...
    'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
set(plot2(4),'DisplayName','1000 Bubbles','Marker','diamond',...
    'Color',[0.611764705882353 0.596078431372549 0.596078431372549]);
% Create ylabel
ylabel('residual log[Err(n)]');

% Create xlabel
xlabel('Number of Iterations [n]');

% grid(axes1,'on');
% Set the remaining axes properties

legend([legend_data1 legend_data2])