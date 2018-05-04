clear all; close all; clc;

if exist('faccept.dat')
    delimiterIn = ' ';
    finaccept = importdata('faccept.dat',delimiterIn);

    delimiterIn = '\t';
    for i=1:2000
        if exist(['test' num2str(i) '.dat'])
            filename = ['test' num2str(i)];
            fintest = importdata([filename '.dat'],delimiterIn);
        end
        
    end

    %%
    delimiterIn = ' ';
    headerlinesIn = 3;
    for i = 1:size(finaccept,1)
        fnameini = strcat('min',num2str(finaccept(i,1)));
        fnamesad = strcat('min',num2str(finaccept(i,2)));
        ini_pos(i) = importdata(fnameini,delimiterIn,headerlinesIn);
        sad_pos(i) = importdata(fnamesad,delimiterIn,headerlinesIn);
    end

    %%
    fid1 = fopen('critwofinal.dat', 'wt'); % Open for writing
    for i = 1:size(finaccept,1)
        fprintf(fid1, '%d %d %d %d %d', fintest.data(i,:), finaccept(i,2));
        fprintf(fid1, '\n');
    end
    fclose(fid1);
end
