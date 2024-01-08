clc;clear all
close all

% load
fin = 'D:\Chronos_Footage\1.1g\output_n';

fps = 3000;
nbins = 100;
dfit = 4;

%% 

flist = dir([fin '\tracer*.mat']);
Nexp = size(flist,1);

for nexp = 1:1
    load([fin '\tracer_STB_' num2str(nexp) '.mat' ]);
    
end
