function LagStatsParticle = LagragianStatParticle(tracklong_particle,Fs,pdf,ndt)

%% Tracer: PDFs of velocity and acceleration 
% Nbins = 256; % bins for PDF
% n = 10; % spacing for PDF

Nbins = pdf.Nbins;
n = pdf.n;

% disp('calculating PDFs')
pdfVx = mkpdf5(tracklong_particle,'Vx',Nbins,n);
pdfVy = mkpdf5(tracklong_particle,'Vy',Nbins,n);
pdfVz = mkpdf5(tracklong_particle,'Vz',Nbins,n);
pdfAx = mkpdf5(tracklong_particle,'Ax',Nbins,2*n);
pdfAy = mkpdf5(tracklong_particle,'Ay',Nbins,2*n);
pdfAz = mkpdf5(tracklong_particle,'Az',Nbins,2*n);

LagStatsParticle.pdfVx = pdfVx;
LagStatsParticle.pdfVy = pdfVy;
LagStatsParticle.pdfVz = pdfVz;
LagStatsParticle.pdfAx = pdfAx;
LagStatsParticle.pdfAy = pdfAy;
LagStatsParticle.pdfAz = pdfAz;

%% Tracer: compute MSD 
% disp('calculating MSDs')
MSDx = structFunc_struct(tracklong_particle,'Xf',2);
MSDy = structFunc_struct(tracklong_particle,'Yf',2);
MSDz = structFunc_struct(tracklong_particle,'Zf',2);

LagStatsParticle.MSDx = MSDx;
LagStatsParticle.MSDy = MSDy;
LagStatsParticle.MSDz = MSDz;
%% Tracer: Lagrangian 2nd SF
% disp('calculating S2L')
[S2Lx, ~, ~]= structFunc_struct(tracklong_particle,'Vx',2);
[S2Ly, ~, ~]= structFunc_struct(tracklong_particle,'Vy',2);
[S2Lz, ~, ~]= structFunc_struct(tracklong_particle,'Vz',2);

LagStatsParticle.S2Lx = S2Lx;
LagStatsParticle.S2Ly = S2Ly;
LagStatsParticle.S2Lz = S2Lz;
%% Tracer: correlation function (fit)

% disp('calculating Correlation functions')
Ruux = xcorr_struct(tracklong_particle,'Vx',1);
Ruuy = xcorr_struct(tracklong_particle,'Vy',1);
Ruuz = xcorr_struct(tracklong_particle,'Vz',1);
Raax = xcorr_struct(tracklong_particle,'Ax',1);
Raay = xcorr_struct(tracklong_particle,'Ay',1);
Raaz = xcorr_struct(tracklong_particle,'Az',1);

LagStatsParticle.Ruux = Ruux;
LagStatsParticle.Ruuy = Ruuy;
LagStatsParticle.Ruuz = Ruuz;
LagStatsParticle.Raax = Raax;
LagStatsParticle.Raay = Raay;
LagStatsParticle.Raaz = Raaz;

%% Tracer: correlation function -- dt method (denosied)
% ndts = [8 8 8]; % start points of correlation function ndt
% ndtl = [10 10 10]; % length of correlation function ndt

ndts = ndt.start;
ndtl = ndt.len;

% disp('calculating Correlation functions -- dt method')
[tau,corrv,corra] = dtCorr(tracklong_particle,ndts,ndtl,Fs);

LagStatsParticle.dtCorrvX.ndts = ndts;
LagStatsParticle.dtCorrvX.ndtl = ndtl;
LagStatsParticle.dtCorrvX.tau = tau.X';
LagStatsParticle.dtCorrvX.corr = corrv.X';

LagStatsParticle.dtCorrvY.ndts = ndts;
LagStatsParticle.dtCorrvY.ndtl = ndtl;
LagStatsParticle.dtCorrvY.tau = tau.Y';
LagStatsParticle.dtCorrvY.corr = corrv.Y';

LagStatsParticle.dtCorrvZ.ndts = ndts;
LagStatsParticle.dtCorrvZ.ndtl = ndtl;
LagStatsParticle.dtCorrvZ.tau = tau.Z';
LagStatsParticle.dtCorrvZ.corr = corrv.Z';

LagStatsParticle.dtCorraX.ndts = ndts;
LagStatsParticle.dtCorraX.ndtl = ndtl;
LagStatsParticle.dtCorraX.tau = tau.X';
LagStatsParticle.dtCorraX.corr = corra.X';

LagStatsParticle.dtCorraY.ndts = ndts;
LagStatsParticle.dtCorraY.ndtl = ndtl;
LagStatsParticle.dtCorraY.tau = tau.Y';
LagStatsParticle.dtCorraY.corr = corra.Y';

LagStatsParticle.dtCorraZ.ndts = ndts;
LagStatsParticle.dtCorraZ.ndtl = ndtl;
LagStatsParticle.dtCorraZ.tau = tau.Z';
LagStatsParticle.dtCorraZ.corr = corra.Z';
clear ndts ndtl