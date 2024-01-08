clear all;clc;close all

%%
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')

%%
Fs=2996;
%%
fname = 'STB';
load(fname);
%%
frame0 = 0;
part(1:frame0-1)=[];
%% Build track structure 
track = part2track(part);

%% Finding long tracks
L = arrayfun(@(X)(numel(X.X)),track);
long_thres = 10;
Ilong = find(L>=long_thres);

%% Find proper filter width
[s(1), m(1), w]=findFilterWidth_PTV(track(Ilong),'X');
[s(2), m(2), w]=findFilterWidth_PTV(track(Ilong),'Y');
[s(3), m(3), w]=findFilterWidth_PTV(track(Ilong),'Z');
%% find optimal filter length
plotLagr_std(w,s,4,mycolormap,0,[])

%% set values for gaussian filter
wopt = 4;
lopt = 20;
%% Estimate smoothed tracks, velocities and accelerations with optimal filter width
tracklong=calcVelLEM(track(Ilong),wopt,lopt,Fs);
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),tracklong)==1);
tracklong = tracklong(Ine);
save('tracks.mat','long_thres','tracklong','wopt','lopt','Fs')

LagFilename = './LagrangianStats.mat';
save(LagFilename,'s','w','wopt','lopt')

%% find optimal ndt for Correlation function -- dt method
nmaxdt = 10;
nmaxtau = 10;
find_optimal_ndt(tracklong,nmaxdt,nmaxtau,'X',1)
find_optimal_ndt(tracklong,nmaxdt,nmaxtau,'Y',1)
find_optimal_ndt(tracklong,nmaxdt,nmaxtau,'Z',1)
clear nmaxtau nmaxtau


%% Lagragian statistics
Nbins = 256; % bins for PDF
n = 10; % spacing for PDF
ndts = [6 6 6]; % start points of correlation function ndt
ndtl = [10 10 10]; % length of correlation function ndt
LagragianStats = LagrangianStats(tracklong,Fs,Nbins,n,ndts,ndtl,1);
clear Nbins n ndts ndtl

%% Lagrangian plots
nCorrFit = [300 20];
if2layers = [1 0];
plotLagr(LagFilename,mycolormap,Fs,nCorrFit,if2layers,1)
clear nCorrFit if2layers

%% 2 point statistics
tic
[eulerStats,pair]= twoPointsEulerianStats_Mica_Speedup(tracklong,[0.5 80],100,0,0);
save('EulerianStats.mat','eulerStats')   
save('EulerianPair.mat','pair','-v7.3')
toc
%% Eulerian plots
plotEuler(eulerStats,Fs,mycolormap,0,1)

%%  mean fields and substract 
dt = [2 4 6 8];
nbins = [20 21 22];
slices.x = [-20 0 20];
slices.y = [0];
slices.z = [-5];
axisrange = [-40 40 -25 30 -30 25];
%%
[gridsV,meanFieldsV,tracklong] = meanFields(tracklong,Fs,dt,nbins,1,1,1);
[gridsVrms,meanFieldsVrms,tracklong] = meanFields(tracklong,Fs,dt,nbins,1,2,1);
[gridsA,meanFieldsA,tracklong] = meanFields(tracklong,Fs,dt,nbins,2,1,1);
[gridsArms,meanFieldsArms,tracklong] = meanFields(tracklong,Fs,dt,nbins,2,2,1);
%%
sliceFields(gridsV,meanFieldsV,slices,axisrange,1,1)
sliceFields(gridsVrms,meanFieldsVrms,slices,axisrange,1,2)
sliceFields(gridsA,meanFieldsA,slices,axisrange,2,1)
sliceFields(gridsArms,meanFieldsArms,slices,axisrange,2,2)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)
% quiverFields(XX,YY,ZZ,mVx,mVy,mVz,axisrange)

%%
clearvars -except tracklong Fs mycolormap

%% 
track_subsMean = trackSubsMean(tracklong);
save('tracks_subsMean.mat','track_subsMean')
clear tracklong
%%
tic
[eulerStats_subsMean,pair_subsMean]= twoPointsEulerianStats_Mica_Speedup(track_subsMean,[0.5 80],100,1,0);
save('EulerianStats_subsMean.mat','eulerStats_subsMean')   
save('EulerianPair_subsMean.mat','pair_subsMean','-v7.3')
toc
%% Eulerian plots
plotEuler(eulerStats_subsMean,Fs,mycolormap,1,1)