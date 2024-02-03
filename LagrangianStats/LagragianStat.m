function LagragianStat(fin,Fs,pdf,ndt,fout,kexp)

load([fin 'tracks_' num2str(kexp) '.mat'],'tracklong_tracers');


% %% Tracer: PDFs of velocity and acceleration 
% % Nbins = 256; % bins for PDF
% % n = 10; % spacing for PDF
% 
% Nbins = pdf.Nbins;
% n = pdf.n;
% 
% % disp('calculating PDFs')
% pdfVx = mkpdf5(tracklong_tracers,'Vx',Nbins,n);
% pdfVy = mkpdf5(tracklong_tracers,'Vy',Nbins,n);
% pdfVz = mkpdf5(tracklong_tracers,'Vz',Nbins,n);
% pdfAx = mkpdf5(tracklong_tracers,'Ax',Nbins,2*n);
% pdfAy = mkpdf5(tracklong_tracers,'Ay',Nbins,2*n);
% pdfAz = mkpdf5(tracklong_tracers,'Az',Nbins,2*n);
% 
% LagrangianStats.pdfVx = pdfVx;
% LagrangianStats.pdfVy = pdfVy;
% LagrangianStats.pdfVz = pdfVz;
% LagrangianStats.pdfAx = pdfAx;
% LagrangianStats.pdfAy = pdfAy;
% LagrangianStats.pdfAz = pdfAz;
% 
% %% Tracer: compute MSD 
% % disp('calculating MSDs')
% MSDx = structFunc_struct(tracklong_tracers,'Xf',2);
% MSDy = structFunc_struct(tracklong_tracers,'Yf',2);
% MSDz = structFunc_struct(tracklong_tracers,'Zf',2);
% 
% LagrangianStats.MSDx = MSDx;
% LagrangianStats.MSDy = MSDy;
% LagrangianStats.MSDz = MSDz;
% %% Tracer: Lagrangian 2nd SF
% % disp('calculating S2L')
% [S2Lx, ~, ~]= structFunc_struct(tracklong_tracers,'Vx',2);
% [S2Ly, ~, ~]= structFunc_struct(tracklong_tracers,'Vy',2);
% [S2Lz, ~, ~]= structFunc_struct(tracklong_tracers,'Vz',2);
% 
% LagrangianStats.S2Lx = S2Lx;
% LagrangianStats.S2Ly = S2Ly;
% LagrangianStats.S2Lz = S2Lz;
% %% Tracer: correlation function (fit)
% 
% % disp('calculating Correlation functions')
% Ruux = xcorr_struct(tracklong_tracers,'Vx',1);
% Ruuy = xcorr_struct(tracklong_tracers,'Vy',1);
% Ruuz = xcorr_struct(tracklong_tracers,'Vz',1);
% Raax = xcorr_struct(tracklong_tracers,'Ax',1);
% Raay = xcorr_struct(tracklong_tracers,'Ay',1);
% Raaz = xcorr_struct(tracklong_tracers,'Az',1);
% 
% LagrangianStats.Ruux = Ruux;
% LagrangianStats.Ruuy = Ruuy;
% LagrangianStats.Ruuz = Ruuz;
% LagrangianStats.Raax = Raax;
% LagrangianStats.Raay = Raay;
% LagrangianStats.Raaz = Raaz;

%% Tracer: correlation function -- dt method (denosied)
% ndts = [8 8 8]; % start points of correlation function ndt
% ndtl = [10 10 10]; % length of correlation function ndt

ndts = ndt.start;
ndtl = ndt.len;

% disp('calculating Correlation functions -- dt method')
[tau,corrv,corra] = dtCorr(tracklong_tracers,ndts,ndtl,Fs);

LagrangianStats.dtCorrvX.ndts = ndts;
LagrangianStats.dtCorrvX.ndtl = ndtl;
LagrangianStats.dtCorrvX.tau = tau.X';
LagrangianStats.dtCorrvX.corr = corrv.X';

LagrangianStats.dtCorrvY.ndts = ndts;
LagrangianStats.dtCorrvY.ndtl = ndtl;
LagrangianStats.dtCorrvY.tau = tau.Y';
LagrangianStats.dtCorrvY.corr = corrv.Y';

LagrangianStats.dtCorrvZ.ndts = ndts;
LagrangianStats.dtCorrvZ.ndtl = ndtl;
LagrangianStats.dtCorrvZ.tau = tau.Z';
LagrangianStats.dtCorrvZ.corr = corrv.Z';

LagrangianStats.dtCorraX.ndts = ndts;
LagrangianStats.dtCorraX.ndtl = ndtl;
LagrangianStats.dtCorraX.tau = tau.X';
LagrangianStats.dtCorraX.corr = corra.X';

LagrangianStats.dtCorraY.ndts = ndts;
LagrangianStats.dtCorraY.ndtl = ndtl;
LagrangianStats.dtCorraY.tau = tau.Y';
LagrangianStats.dtCorraY.corr = corra.Y';

LagrangianStats.dtCorraZ.ndts = ndts;
LagrangianStats.dtCorraZ.ndtl = ndtl;
LagrangianStats.dtCorraZ.tau = tau.Z';
LagrangianStats.dtCorraZ.corr = corra.Z';
clear ndts ndtl


%% Tracer: save Lagrangian Stats
% save([fout 'LagrangianStats_' num2str(kexp) '.mat'],'LagrangianStats','-append')
save([fout 'Corr_' num2str(kexp) '.mat'],'LagrangianStats')
% disp('---- Lagrangian Stats done ----')