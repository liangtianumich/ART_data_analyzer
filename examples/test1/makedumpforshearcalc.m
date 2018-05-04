clear all; close all; clc;

delimiterIn = ' ';

if exist('acceptlist.dat')
faccept = importdata('acceptlist.dat',delimiterIn);

for ii = 1:size(faccept)
%% for minimum position
% for j = 1000:1009 
%     if j == faccept(ii,1)
%         filename = ['min' num2str(j)];
%         refilename = [num2str(ii) '-' filename];
% %-----Read data------------------------------------------------------------
% d = 2000; % number of atoms  
% [ele, data1, data2, data3] = textread(filename,'%n%f%f%f', 'headerlines', 3);
% for i=1:1:d
%     no(i) = i;
% end
% M = [no', ele, data1, data2, data3];
% %-----write data-----------------------------------------------------------
% fileOut = ['M' refilename '.dump'];
% fid=fopen(fileOut,'wt'); [m,n]=size(M);
% fprintf(fid,'%s\n', 'ITEM: TIMESTEP');
% fprintf(fid,'%s\n', '0');
% fprintf(fid,'%s\n', 'ITEM: NUMBER OF ATOMS');
% fprintf(fid,'%d\n', d);
% fprintf(fid,'%s\n', 'ITEM: BOX BOUNDS pp pp pp');
% fprintf(fid,'%f\t', 2.9987472714648433e-01); fprintf(fid,'%f\n', 3.2130125270428323e+01);
% fprintf(fid,'%f\t', 2.9987472714648433e-01); fprintf(fid,'%f\n', 3.2130125270428323e+01);
% fprintf(fid,'%f\t', 2.9987472714648433e-01); fprintf(fid,'%f\n', 3.2130125270428323e+01);
% fprintf(fid,'%s\n', 'ITEM: ATOMS id type x y z');
% for i=1:1:m
%      fprintf(fid,'%d\t',M(i,1));
%      fprintf(fid,'%d\t',M(i,2));
%      fprintf(fid,'%8.6f\t',M(i,3));
%      fprintf(fid,'%8.6f\t',M(i,4));
%      fprintf(fid,'%8.6f\n',M(i,5));
% end
% end
% end

%% for saddle position
% for k = 1001:1010 
%     if k == faccept(ii,2)
%         filename1 = ['sad' num2str(k)];
%         refilename1 = [num2str(ii) '-' filename1];
% 
% %-----Read data------------------------------------------------------------
% d = 2000; % number of atoms  
% [ele, data1, data2, data3] = textread(filename1,'%n%f%f%f', 'headerlines', 3);
% for i=1:1:d
%     no(i) = i;
% end
% M = [no', ele, data1, data2, data3];
% 
% %-----write data-----------------------------------------------------------
% fileOut = ['M' refilename1 '.dump'];
% fid1=fopen(fileOut,'wt'); [m,n]=size(M);
% fprintf(fid1,'%s\n', 'ITEM: TIMESTEP');
% fprintf(fid1,'%s\n', '0');
% fprintf(fid1,'%s\n', 'ITEM: NUMBER OF ATOMS');
% fprintf(fid1,'%d\n', d);
% fprintf(fid1,'%s\n', 'ITEM: BOX BOUNDS pp pp pp');
% fprintf(fid1,'%f\t', 2.9987472714648433e-01); fprintf(fid1,'%f\n', 3.2130125270428323e+01);
% fprintf(fid1,'%f\t', 2.9987472714648433e-01); fprintf(fid1,'%f\n', 3.2130125270428323e+01);
% fprintf(fid1,'%f\t', 2.9987472714648433e-01); fprintf(fid1,'%f\n', 3.2130125270428323e+01);
% fprintf(fid1,'%s\n', 'ITEM: ATOMS id type x y z');
% for i=1:1:m
%      fprintf(fid1,'%d\t',M(i,1));
%      fprintf(fid1,'%d\t',M(i,2));
%      fprintf(fid1,'%8.6f\t',M(i,3));
%      fprintf(fid1,'%8.6f\t',M(i,4));
%      fprintf(fid1,'%8.6f\n',M(i,5));
% end
% end
% end

%% for final position
for k = 1001:1010 
   if k == faccept(ii,2)
       filename2 = ['min' num2str(k)];
       refilename2 = [num2str(ii) '-' filename2];

%-----Read data------------------------------------------------------------
d = 2000; % number of atoms  
[ele, data1, data2, data3] = textread(filename2,'%n%f%f%f', 'headerlines', 3);
for i=1:1:d
   no(i) = i;
end
M = [no', ele, data1, data2, data3];

%-----write data-----------------------------------------------------------
fileOut = ['F' refilename2 '.dump'];
fid2=fopen(fileOut,'wt'); [m,n]=size(M);
fprintf(fid2,'%s\n', 'ITEM: TIMESTEP');
fprintf(fid2,'%s\n', '0');
fprintf(fid2,'%s\n', 'ITEM: NUMBER OF ATOMS');
fprintf(fid2,'%d\n', d);
fprintf(fid2,'%s\n', 'ITEM: BOX BOUNDS pp pp pp');
fprintf(fid2,'%f\t', 2.9987472714648433e-01); fprintf(fid2,'%f\n', 3.2130125270428323e+01);
fprintf(fid2,'%f\t', 2.9987472714648433e-01); fprintf(fid2,'%f\n', 3.2130125270428323e+01);
fprintf(fid2,'%f\t', 2.9987472714648433e-01); fprintf(fid2,'%f\n', 3.2130125270428323e+01);
fprintf(fid2,'%s\n', 'ITEM: ATOMS id type x y z');
for i=1:1:m
    fprintf(fid2,'%d\t',M(i,1));
    fprintf(fid2,'%d\t',M(i,2));
    fprintf(fid2,'%8.6f\t',M(i,3));
    fprintf(fid2,'%8.6f\t',M(i,4));
    fprintf(fid2,'%8.6f\n',M(i,5));
end
end
end



end
end
