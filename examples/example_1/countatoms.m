clear all; close all; clc;

sheardata = importdata('teststrain.dat'); %('Cu64Zr36_10_11_shear.dat');  %read data from files
minpos = importdata('testmin.dat');
sadpos = importdata('testsad.dat');

for i=1:length(minpos)
    dist(i,1) = sqrt((sadpos(i,3)-minpos(i,3))^(2)+(sadpos(i,4)-minpos(i,4))^(2)+ ...
        (sadpos(i,5)-minpos(i,5))^(2));
end


%%
X = dist(:,1);
Y = sheardata(:,7);

% figure(1)
% scatter(X,Y)
% xlabel('displacement')
% ylabel('shear strain')
% hold on;

% linear fitting
 f = polyfit(X,Y,1);
 y_fit = polyval(f,X);
% plot(X,y_fit,'r')

%%
%rd = 7.6*std(Y);
%y1 = y_fit + rd;
%y2 = y_fit - rd;
 
% figure(2)
% y3 = y_fit+rd;
% y4 = y_fit-rd;
% scatter(X,Y)
% hold on
% plot(X,y_fit,'r')
% plot(X,y3)
% plot(X,y4) 

% count = 0;
% for i = 1:length(sheardata(:,7))
%     if (Y(i) > y1(i) || Y(i) < y2(i))
%         count = count+1;
%     end
% end

%%
P = [dist(:,1) sheardata(:,7)];
a = f(1,1); % 0.0909020588381900;
b = f(1,2); %0.00295249493742894;

for i=1:length(dist(:,1))
    d(i,1) = abs(a*P(i,1) + b*P(i,2))/sqrt(a^2 + b^2);
end

deviation = [0.2; 0.4; 0.6; 0.8; 1; 1.2; 1.4; 1.6; 1.8; 2];
numatoms1 = 0; numatoms2 = 0; numatoms3 = 0; numatoms4 = 0; numatoms5 = 0;
numatoms6 = 0; numatoms7 = 0; numatoms8 = 0; numatoms9 = 0; numatoms10 =0;
for i=1:length(d(:,1))
    if (0 < d(i,1)) && (d(i,1) < deviation(1))
        numatoms1 = numatoms1 + 1;
        count1(i,1) = numatoms1;
        count1(i,2) = i;
    elseif (deviation(1) < d(i,1)) && (d(i,1) < deviation(2))
        numatoms2 = numatoms2 + 1;
        count2(i,1) = numatoms2;
        count2(i,2) = i;
    elseif (deviation(2) < d(i,1)) && (d(i,1) < deviation(3))
        numatoms3 = numatoms3 + 1;
        count3(i,1) = numatoms3;
        count3(i,2) = i;
    elseif (deviation(3) < d(i,1)) && (d(i,1) < deviation(4))
        numatoms4 = numatoms4 + 1;
        count4(i,1) = numatoms4;
        count4(i,2) = i;
    elseif (deviation(4) < d(i,1)) && (d(i,1) < deviation(5))
        numatoms5 = numatoms5 + 1;
        count5(i,1) = numatoms5;
        count5(i,2) = i;
    elseif (deviation(5) < d(i,1)) && (d(i,1) < deviation(6))
        numatoms6 = numatoms6 + 1;
        count6(i,1) = numatoms6;
        count6(i,2) = i;
    elseif (deviation(6) < d(i,1)) && (d(i,1) < deviation(7))
        numatoms7 = numatoms7 + 1;
        count7(i,1) = numatoms7;
        count7(i,2) = i;
    elseif (deviation(7) < d(i,1)) && (d(i,1) < deviation(8))
        numatoms8 = numatoms8 + 1;
        count8(i,1) = numatoms8;
        count8(i,2) = i;
    elseif (deviation(8) < d(i,1)) && (d(i,1) < deviation(9))
        numatoms9 = numatoms9 + 1;
        count9(i,1) = numatoms9;
        count9(i,2) = i;
    else
        numatoms10 = numatoms10 + 1;
        count10(i,1) = numatoms10;
        count10(i,2) = i;
    end
end

atomcount = [numatoms1; numatoms2; numatoms3; numatoms4; numatoms5; ...
    numatoms6; numatoms7; numatoms8; numatoms9; numatoms10];

if (numatoms1 ~= 0)
    count1 = count1(any(count1,2),:);
end
if (numatoms2 ~= 0)
    count2 = count2(any(count2,2),:);
end
if (numatoms3 ~= 0)
   count3 = count3(any(count3,2),:);
end
if (numatoms4 ~= 0)
    count4 = count4(any(count4,2),:);
end
if (numatoms5 ~= 0)
    count5 = count5(any(count5,2),:);
end
if (numatoms6 ~= 0)
    count6 = count6(any(count6,2),:);
end
if (numatoms7 ~= 0)
    count7 = count7(any(count7,2),:);
end
if (numatoms8 ~= 0)
    count8 = count8(any(count8,2),:);
end
if (numatoms9 ~= 0)
    count9 = count9(any(count9,2),:);
end
if (numatoms10 ~= 0)
    count10 = count10(any(count10,2),:);
end

triger = numatoms5 + numatoms6 + numatoms7 + numatoms8 + ...
    numatoms9 + numatoms10; 

% figure(3)
% scatter(deviation, atomcount)
% xlim([0 2])
% ylim([1 2000])

%%
fid = fopen('numatominvolved.dat', 'wt'); % Open for writing
fprintf(fid, '%d', triger);
fclose(fid);

%if (numatoms3~=0)
%    fid = fopen('trigeratomid3.dat','wt');
%    fprintf(fid, '%d \n', count3(:,2)); 
%    fclose(fid);
%end

%if (numatoms4~=0)
%    fid1 = fopen('trigeratomid4.dat','wt');
%    fprintf(fid1, '%d \n', count4(:,2)); 
%    fclose(fid1);
%end

%if (numatoms5~=0)
%    fid2 = fopen('trigeratomid5.dat','wt');
%    fprintf(fid2, '%d \n', count5(:,2)); 
%    fclose(fid2);
%end

%if (numatoms6~=0)
%    fid3 = fopen('trigeratomid6.dat','wt');
%    fprintf(fid3, '%d \n', count6(:,2)); 
%    fclose(fid3);
%end

%if (numatoms7~=0)
%    fid4 = fopen('trigeratomid7.dat','wt');
%    fprintf(fid4, '%d \n', count7(:,2)); 
%    fclose(fid4);
%end

%if (numatoms8~=0)
%    fid5 = fopen('trigeratomid8.dat','wt');
%    fprintf(fid5, '%d \n', count8(:,2)); 
%    fclose(fid5);
%end

%if (numatoms9~=0)
%    fid6 = fopen('trigeratomid9.dat','wt');
%    fprintf(fid6, '%d \n', count9(:,2)); 
%    fclose(fid6);
%end

%if (numatoms10~=0)
%    fid7 = fopen('trigeratomid10.dat','wt');
%    fprintf(fid7, '%d \n', count10(:,2)); 
%    fclose(fid7);
%end
