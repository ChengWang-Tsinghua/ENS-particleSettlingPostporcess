function [averLS,LS] = averLagrangianStats(lagrangianStats)

% LS.pdfnVx = cell2mat(arrayfun(@(X)(X.pdfn),lagrangianStats.pdfVx,'UniformOutput',false));
% LS.pdfnVy = cell2mat(arrayfun(@(X)(X.pdfn),lagrangianStats.pdfVy,'UniformOutput',false));
% LS.pdfnVz = cell2mat(arrayfun(@(X)(X.pdfn),lagrangianStats.pdfVz,'UniformOutput',false));
% LS.xpdfnVx = cell2mat(arrayfun(@(X)(X.xpdfn),lagrangianStats.pdfVx,'UniformOutput',false));
% LS.xpdfnVy = cell2mat(arrayfun(@(X)(X.xpdfn),lagrangianStats.pdfVy,'UniformOutput',false));
% LS.xpdfnVz = cell2mat(arrayfun(@(X)(X.xpdfn),lagrangianStats.pdfVz,'UniformOutput',false));
% 
% LS.pdfnAx = cell2mat(arrayfun(@(X)(X.pdfn),lagrangianStats.pdfAx,'UniformOutput',false));
% LS.pdfnAy = cell2mat(arrayfun(@(X)(X.pdfn),lagrangianStats.pdfAy,'UniformOutput',false));
% LS.pdfnAz = cell2mat(arrayfun(@(X)(X.pdfn),lagrangianStats.pdfAz,'UniformOutput',false));
% LS.xpdfnAx = cell2mat(arrayfun(@(X)(X.xpdfn),lagrangianStats.pdfAx,'UniformOutput',false));
% LS.xpdfnAy = cell2mat(arrayfun(@(X)(X.xpdfn),lagrangianStats.pdfAy,'UniformOutput',false));
% LS.xpdfnAz = cell2mat(arrayfun(@(X)(X.xpdfn),lagrangianStats.pdfAz,'UniformOutput',false));
% 
% LS.tauMSDx = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.MSDx,'UniformOutput',false));
% LS.tauMSDy = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.MSDy,'UniformOutput',false));
% LS.tauMSDz = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.MSDz,'UniformOutput',false));
% LS.meanMSDx = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.MSDx,'UniformOutput',false));
% LS.meanMSDy = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.MSDy,'UniformOutput',false));
% LS.meanMSDz = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.MSDz,'UniformOutput',false));
% 
% LS.tauS2Lx = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.S2Lx,'UniformOutput',false));
% LS.tauS2Ly = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.S2Ly,'UniformOutput',false));
% LS.tauS2Lz = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.S2Lz,'UniformOutput',false));
% LS.meanS2Lx = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.S2Lx,'UniformOutput',false));
% LS.meanS2Ly = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.S2Ly,'UniformOutput',false));
% LS.meanS2Lz = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.S2Lz,'UniformOutput',false));
% 
% LS.tauRuux = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.Ruux,'UniformOutput',false));
% LS.tauRuuy = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.Ruuy,'UniformOutput',false));
% LS.tauRuuz = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.Ruuz,'UniformOutput',false));
% LS.meanRuux = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.Ruux,'UniformOutput',false));
% LS.meanRuuy = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.Ruuy,'UniformOutput',false));
% LS.meanRuuz = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.Ruuz,'UniformOutput',false));
% 
% LS.tauRaax = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.Raax,'UniformOutput',false));
% LS.tauRaay = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.Raay,'UniformOutput',false));
% LS.tauRaaz = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.Raaz,'UniformOutput',false));
% LS.meanRaax = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.Raax,'UniformOutput',false));
% LS.meanRaay = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.Raay,'UniformOutput',false));
% LS.meanRaaz = vertcat_pad(arrayfun(@(X)(X.mean),lagrangianStats.Raaz,'UniformOutput',false));

LS.taudtCorrvx = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.dtCorrvX,'UniformOutput',false));
LS.taudtCorrvy = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.dtCorrvY,'UniformOutput',false));
LS.taudtCorrvz = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.dtCorrvZ,'UniformOutput',false));
LS.corrdtCorrvx = vertcat_pad(arrayfun(@(X)(X.corr),lagrangianStats.dtCorrvX,'UniformOutput',false));
LS.corrdtCorrvy = vertcat_pad(arrayfun(@(X)(X.corr),lagrangianStats.dtCorrvY,'UniformOutput',false));
LS.corrdtCorrvz = vertcat_pad(arrayfun(@(X)(X.corr),lagrangianStats.dtCorrvZ,'UniformOutput',false));

LS.taudtCorrax = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.dtCorraX,'UniformOutput',false));
LS.taudtCorray = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.dtCorraY,'UniformOutput',false));
LS.taudtCorraz = vertcat_pad(arrayfun(@(X)(X.tau),lagrangianStats.dtCorraZ,'UniformOutput',false));
LS.corrdtCorrax = vertcat_pad(arrayfun(@(X)(X.corr),lagrangianStats.dtCorraX,'UniformOutput',false));
LS.corrdtCorray = vertcat_pad(arrayfun(@(X)(X.corr),lagrangianStats.dtCorraY,'UniformOutput',false));
LS.corrdtCorraz = vertcat_pad(arrayfun(@(X)(X.corr),lagrangianStats.dtCorraZ,'UniformOutput',false));

fieldname = fieldnames(LS);
for i = 1:numel(fieldname)
    averLS.(fieldname{i}) = mean(LS.(fieldname{i}),1,"omitnan");
end