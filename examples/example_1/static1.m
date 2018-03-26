clear all; close all; clc;

delimiterIn = '\t';
for i=1:2000
    if exist(['atominfo' num2str(i) '.dat'])
       filename = ['atominfo' num2str(i) '.dat'];
       data = importdata([filename ],delimiterIn);
    end
end
num = length(data); 
%%
atomdata = zeros();
for i = 1:num-1
    if data(i+1,2)-data(i,2) == 1
      atomdata(i,1) = data(i,1);
      atomdata(i,2) = data(i,2);
    end
end
atomdata(num,1) = data(num,1);
atomdata(num,2) = data(num,2);

atomdata = atomdata(any(atomdata,2),:); % remove zeros from matrix

for i = 1:length(atomdata)
    %if atomdata(i,2) == i
    atomdata(i,3) = i+1000;
end
        

%%
fid = fopen('atomdata.dat', 'wt');
for i = 1:length(atomdata)
    fprintf(fid, '%d %d %d %s', atomdata(i,:), filename);
    fprintf(fid, '\n');
end
fclose(fid)