clear all; close all; clc;

if exist('faccept.dat')
    delimiterIn = '\t';
    delimiterIn1 = ' ';
    ini_pos = importdata('testmin.dat',delimiterIn1);
    sad_pos = importdata('testsad.dat',delimiterIn);
    av_atoms = importdata('testav.dat');
    
    %%
    for i = 1:length(ini_pos)
        dist(i,1) = (sad_pos(i,3)-ini_pos(i,3))^2 + ...
       (sad_pos(i,4)-ini_pos(i,4))^2 + (sad_pos(i,5)-ini_pos(i,5))^2;
    end
    
    msd = sum(dist)/length(dist);
    %ave_av = sum(av_atoms(:,3))/length(av_atoms);
    
    %a = (ave_av)^(1/3);
    for i = 1:length(av_atoms)
        a(i,1) = (av_atoms(i,3))^(1/3);
    end
        
    %vflex = msd/(a^2)*ave_av;
   for i=1:length(dist)
       vflex(i,1) = msd/(a(i,1)^2)*av_atoms(i,3);
   end
   
   ave_vflex = sum(vflex)/length(vflex);
   
   %%
   for i=1:length(ini_pos)
       if (vflex(i,1) > max(vflex)*0.97)
           topvflex(i,1) = vflex(i,1);
       elseif (vflex(i,1) < min(vflex)*1.03)
           botvflex(i,1) = vflex(i,1);
       end
   end
   
   topvflex = topvflex(any(topvflex,2),:);
   botvflex = botvflex(any(botvflex,2),:);
   ave_topvflex = mean(topvflex);
   ave_botvflex = mean(botvflex);
  %% 
  C = 9/(4*(pi^2))*((4*pi)/9)^(2/3)*(17/8)^(2/3);
  kb = 1.38064852*10^(-23); % m2 kg s-2 K-1
  T = 300; % [K]
  
  for i=1:length(vflex)
      G(i,1) = C*kb*T/(vflex(i,1)*(10^(-10)^2));
  end
  
  for i=1:length(topvflex)
      topG(i,1) = C*kb*T/(topvflex(i,1)*(10^(-10)^2));
  end
  
  for i=1:length(botvflex)
      botG(i,1) = C*kb*T/(botvflex(i,1)*(10^(-10)^2));
  end
  ave_G = sum(G)/length(G);
  ave_topG = mean(topG);
  ave_botG = mean(botG);
  
    %%
%     fid = fopen('vflex.dat', 'wt'); % Open for writing
%     fprintf(fid, '%f ', vflex);
%     fprintf(fid, '\n');
%     fclose(fid);

    fid = fopen('vflex.dat', 'wt'); % Open for writing
    for i=1:length(vflex)
        fprintf(fid, '%f ', vflex(i,1));
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    fid1 = fopen('vflexaverage.dat', 'wt'); % Open for writing
    fprintf(fid1, '%f ', ave_vflex);
    fclose(fid1);
    
    fid2 = fopen('G.dat', 'wt'); % Open for writing
    for i=1:length(G)
        fprintf(fid2, '%f ', G(i,1));
        fprintf(fid2, '\n');
    end
    fclose(fid2);
    
    fid3 = fopen('Gaverage.dat', 'wt'); % Open for writing
    fprintf(fid3, '%f ', ave_G);
    fclose(fid3);
    
    fid4 = fopen('topvflexaverage.dat', 'wt'); % Open for writing
    fprintf(fid4, '%f ', ave_topvflex);
    fclose(fid4);
    
    fid5 = fopen('botvflexaverage.dat', 'wt'); % Open for writing
    fprintf(fid5, '%f ', ave_botvflex);
    fclose(fid5);
end

%%
% figure(1)
% [N1, X1] = hist(vflex(:,1), 50);
% y1D = N1./sum(N1);
% h1D = bar(X1, y1D, 1);
% ylabel('Frequency')
% xlabel('V_{flex}')
% 
% figure(2)
% [N1, X1] = hist(G(:,1), 50);
% y1D = N1./sum(N1);
% h1D = bar(X1, y1D, 1);
% ylabel('Frequency')
% xlabel('Shear modulus')
% 
% figure(3)
% plot(vflex,G)
% 
% figure(4)
% [N1, X1] = hist(topvflex(:,1), 50);
% y1D = N1./sum(N1);
% h1D = bar(X1, y1D, 1);
% ylabel('Frequency')
% xlabel('V_{flex}')
% 
% figure(5)
% [N1, X1] = hist(botvflex(:,1), 50);
% y1D = N1./sum(N1);
% h1D = bar(X1, y1D, 1);
% ylabel('Frequency')
% xlabel('V_{flex}')
