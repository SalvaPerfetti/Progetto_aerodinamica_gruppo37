%% TRASFORMO FILE DA TXT A TABELLE MATLAB
clc
clear
close all
savefile= 'NACA2412.txt';
abc=fopen(savefile);
data=textscan(abc, '%f %f','CollectOutput',1,...
    'Delimiter','','HeaderLines',0);
fclose(abc);
x=data{1}(:,1);
y=data{1}(:,2);