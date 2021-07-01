clear all 
clc

f0                  = 1.0;                                                  % fundamental frequency
f0                  = f0*1000000; % freq in Hz
c0                  = 1482;
% Bubble_In = 1;
filename = '1\BUBBLES\output_file_bubble.bin';
file_str = split(filename,'\');
% fileID = fopen(filename,'r');
% formatSpec = '%f %f %f %f %f %f';
% sizeA = [6 Inf];
% A = fscanf(fileID,formatSpec,sizeA);
BubbleN = str2double(file_str(1));
% iDimT = size(A,2)\BubbleN;
% A = A(:,1+(Bubble_In-1)*iDimT:Bubble_In*iDimT);
% time = A(1,:)';
% pressure = A(2,:)';
% R =A(3,:)';
% R_d = A(4,:)';
% R_dd = A(5,:)';
% V_dd_norm = A(6,:)';
% plot(R)


for i =0:53
if i<=9 filename2 = [int2str(BubbleN),'\BUBBLES\output_file_00',int2str(i),'.bin']; else filename2 = [int2str(BubbleN),'\BUBBLES\output_file_0',int2str(i),'.bin']; end
% filename2 = ['H:\Desktop\INCS3D_Agis_v4_20200212_many_bubble\Source_Code\Project_file\output_file_00',int2str(i),'.bin'];
file_str = split(filename2,'\');
fileID2 = fopen(filename2,'r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
ALoc2= fscanf(fileID2,formatSpec,sizeA);
Aloc2N  = ALoc2*c0*1e+3\f0;
% BubbleN = str2double(file_str(1));
% Bubble_In =1;
% Loc(i+1,:) = ALoc(Bubble_In,:);
end


for i =0:53
if i<=9 filename3 = [int2str(BubbleN),'\BUBBLES\output_file_contrast_00',int2str(i),'.bin']; else filename3 = [int2str(BubbleN),'\BUBBLES\output_file_contrast_0',int2str(i),'.bin']; end
% filename3 = ['H:\Desktop\INCS3D_Agis_v4_20200212_many_bubble\Source_Code\Project_file\output_file_contrast_00',int2str(i),'.bin'];
file_str = split(filename3,'\');
fileID3 = fopen(filename3,'r');
formatSpec = '%f';
sizeA = [1 Inf];
ACon2 = fscanf(fileID3,formatSpec,sizeA);
% BubbleN = str2double(file_str(1));
% Bubble_In =1;
% Con(i+1,:) = ALoc(Bubble_In*iDimT:(Bubble_In+1)*iDimT);
end