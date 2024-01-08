clear;clc;close all

%%
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
%% exp_parameters
Fs=4230;
mp = 0.0048e-3;%kg
dp = 1.0094e-3;%m
rhof = 1000; %kg/m3
nv = 8.532e-7;%m2/s, @ 27 celsius
geff = 1; %ratio to 9.8m/s2

%%
fname = 'STB';
load(fname);

%% remove Lavision velocity&acce. fields
fields = {'Vx','Vy','Vz','V','Ax','Ay','Az','A'};
part = rmfield(part,fields);
clear fields fname
%% Build track structure 
track = part2track(part);

%% Finding long tracks
L = arrayfun(@(X)(numel(X.X)),track);
long_thres = 10;
Ilong = find(L>=long_thres);
clear L
%% Find proper filter width
[s(1), m(1), w]=findFilterWidth_PTV(track(Ilong),'X');
[s(2), m(2), w]=findFilterWidth_PTV(track(Ilong),'Y');
[s(3), m(3), w]=findFilterWidth_PTV(track(Ilong),'Z');
%% find optimal filter length
plotLagr_std(w,s,12,mycolormap,0)

%% set values for gaussian filter
wopt = 12;
lopt = 30;
%% Estimate smoothed tracks, velocities and accelerations with optimal filter width
tracklong=calcVelLEM(track(Ilong),wopt,lopt,Fs);
Ine=arrayfun(@(X)(~isempty(X.Vx)),tracklong)==1;
tracklong = tracklong(Ine);
save('tracks.mat','long_thres','tracklong','wopt','lopt','Fs')

LagFilename = './LagrangianStats.mat';
save(LagFilename,'s','w','wopt','lopt')

clear long_thres track Ilong s m w wopt lopt Ine part

%% traj
part = track2part(tracklong,{'Tf','Xf','Yf','Zf','Vx','Vy','Vz','Ax','Ay','Az'},1);

%% plot traj
fout_traj = './Figures_traj';
titlestring = '$1g-without Turb.$';
axisrange_traj = [-5 25 -18 10 -40 40];
axisrange_vat = [0.05 0.3];
plot_traj(part,titlestring,axisrange_traj,1,fout_traj)
plot_vadiffusion(part,axisrange_vat,Fs,1,fout_traj)
clear fout_traj titlestring axisrange_traj axisrange_vat

%%
taup = 0.10;
params = exp_param(Fs,mp,dp,rhof,nv,geff,taup)
save('./expParams.mat','params');
clear mp dp rhof nv geff taup params
%% find optimal ndt for Correlation function -- dt method
nmaxdt = 30;
nmaxtau = 10;
find_optimal_ndt(tracklong,nmaxdt,nmaxtau,'X',1)
find_optimal_ndt(tracklong,nmaxdt,nmaxtau,'Y',1)
find_optimal_ndt(tracklong,nmaxdt,nmaxtau,'Z',1)
clear nmaxtau nmaxtau

%% Lagragian statistics
Nbins = 256; % bins for PDF
n = 10; % spacing for PDF
ndts = [20 20 20]; % start points of correlation function ndt
ndtl = [10 10 10]; % length of correlation function ndt
LagragianStats = LagrangianStats(tracklong,Fs,Nbins,n,ndts,ndtl,1);
clear Nbins n ndts ndtl

%% Lagrangian plots
close all
nCorrFit = [300 20];
if2layers = [1 0];
plotLagr(LagFilename,mycolormap,Fs,nCorrFit,if2layers,1)
clear nCorrFit if2layers
