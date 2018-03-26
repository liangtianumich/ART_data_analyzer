clear all; close all; clc;

if exist('acceptlist.dat')
acceptinfo = importdata('acceptlist.dat');

fileID = fopen('atomdata.dat');
atominfo = textscan(fileID,'%d %d %d %s');
fclose(fileID);

for i=1:length(acceptinfo(:,1))
    for j=1:length(atominfo{1,3})
        if (acceptinfo(i,2)==atominfo{1,3}(j))
            dispatom(i,1) = atominfo{1,1}(j);
        end
    end
end

%%
 fid = fopen('displacedatom.dat', 'wt'); % Open for writing
 for i=1:length(dispatom)
     fprintf(fid, '%d', dispatom(i,1));
     fprintf(fid, '\n');
 end
 fclose(fid);
end
