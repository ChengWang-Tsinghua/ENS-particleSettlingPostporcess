clc;clear all
close all

%%
addpath(genpath('C:\Users\Gsu\Desktop\SDT_EXP'));
mkdir('./Figures')

%% exp_parameters
Fs=2996;
mp = 0.0048e-3;%kg
dp = 1.0094e-3;%m
rhof = 1000; %kg/m3
nv = 8.532e-7;%m2/s, @ 27 celsius
geff = 1; %ratio to 9.8m/s2
params = exp_param(Fs,mp,dp,rhof,nv,geff)
save('./expParams.mat','params');
%%
load('.\STB.mat')

nbins = 256;
dfit = 3;

%%
tmin = 0.07;
tmax = 0.3;
tfmax = 0.14;
tinterp = tmin:0.0001:tmax;
tfinterp = tmin:0.0001:tfmax;
%%
Nframemax = 1e6;
NTrackmax = 1e9;
nexp= fix(part.T/Nframemax)+1;
t = mod(part.T,Nframemax)/fps;
uni_exp = unique(nexp);
Nexp = numel(uni_exp);


%%
mycolormap = mycolor('#063970','#eeeee4','#e28743');%('#063970','#eeeee4','#e28743')
mycolorind = 1:floor(size(mycolormap,1)/Nexp):size(mycolormap,1);
patchcolor1 = mycolormap(size(mycolormap,1)/2,:);
patchcolor2 = mycolor('#a155b9');
patchalpha1 = 0.4;
patchalpha2 = 0.3;

%% lavision output Velocity
Vinterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% vx
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Vinterp.x = [Vinterp.x; interp1(t(ind),part.Vx(ind),tinterp)];
    plot(t(ind),part.Vx(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.x,"omitnan") + std(Vinterp.x,"omitnan");
lowbound = mean(Vinterp.x,"omitnan") - std(Vinterp.x,"omitnan");
hm = plot(tinterp,mean(Vinterp.x,"omitnan"),'k--',LineWidth=2);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-0.35 0.35])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
legend([hm,hs],'$\langle V_i \rangle$','$\langle V_i \rangle \pm \sigma (V_i)$','interpreter','latex',Location='south',FontSize=12);
title('$Velocity-Lavision$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$V_x (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% vy
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Vinterp.y = [Vinterp.y; interp1(t(ind),part.Vy(ind),tinterp)];
    plot(t(ind),part.Vy(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.y,"omitnan") + std(Vinterp.y,"omitnan");
lowbound = mean(Vinterp.y,"omitnan") - std(Vinterp.y,"omitnan");
hm = plot(tinterp,mean(Vinterp.y,"omitnan"),'k--',LineWidth=2);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-0.7 0])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$V_y (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% vz
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Vinterp.z = [Vinterp.z; interp1(t(ind),part.Vz(ind),tinterp)];
    plot(t(ind),part.Vz(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.z,"omitnan") + std(Vinterp.z,"omitnan");
lowbound = mean(Vinterp.z,"omitnan") - std(Vinterp.z,"omitnan");
hm = plot(tinterp,mean(Vinterp.z,"omitnan"),'k--',LineWidth=2);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-0.35 0.35])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$V_z (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim([-inf,0.35])

%
linkaxes(tf.Children,'x')
colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.')
set(col,'position',[0.42 0.34 0.2 0.02])

savefig('./Figures/V')
saveas(gcf,'./Figures/V','png')
saveas(gcf,'./Figures/V','pdf')

%% lavision output Acceleration
Ainterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% ax
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Ainterp.x = [Ainterp.x; interp1(t(ind),part.Ax(ind),tinterp)];
    plot(t(ind),part.Ax(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.x,"omitnan") + std(Ainterp.x,"omitnan");
lowbound = mean(Ainterp.x,"omitnan") - std(Ainterp.x,"omitnan");
hm = plot(tinterp,mean(Ainterp.x,"omitnan"),'k--',LineWidth=2);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-250 250])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
legend([hm,hs],'$\langle A_i \rangle$','$\langle A_i \rangle \pm \sigma (A_i)$','interpreter','latex',Location='south',FontSize=12);
title('$Acce.-Lavision$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$A_x (m \cdot s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% ay
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Ainterp.y = [Ainterp.y; interp1(t(ind),part.Ay(ind),tinterp)];
    plot(t(ind),part.Ax(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.y,"omitnan") + std(Ainterp.y,"omitnan");
lowbound = mean(Ainterp.y,"omitnan") - std(Ainterp.y,"omitnan");
hm = plot(tinterp,mean(Ainterp.y,"omitnan"),'k--',LineWidth=2);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-250 250])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$A_y (m \cdot s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% az
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Ainterp.z = [Ainterp.z; interp1(t(ind),part.Az(ind),tinterp)];
    plot(t(ind),part.Ax(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.z,"omitnan") + std(Ainterp.z,"omitnan");
lowbound = mean(Ainterp.z,"omitnan") - std(Ainterp.z,"omitnan");
hm = plot(tinterp,mean(Ainterp.z,"omitnan"),'k--',LineWidth=2);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-250 250])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$A_z (m \cdot s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim([-inf,0.35])

%
linkaxes(tf.Children,'x')
colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.');
set(col,'position',[0.42 0.34 0.2 0.02])

savefig('./Figures/A')
saveas(gcf,'./Figures/A','png')
saveas(gcf,'./Figures/A','pdf')

%% fit 
fitpart = struct('T',[],'x',[],'y',[],'z',[],'Vx',[],'Vy',[],'Vz',[],'Ax',[],'Ay',[],'Az',[]);
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    fitpart.t(ind) = t(ind);
    
    % fit x,y,z using polynomial
%     pfit.x= polyfit(fitpart.t(ind),fitpart.X(ind),dfit);
%     pfit.y= polyfit(fitpart.t(ind),fitpart.Y(ind),dfit);
%     pfit.z= polyfit(fitpart.t(ind),fitpart.Z(ind),dfit);
% 
%     fitpart.X(ind) = polyval(pfit.x,fitpart.t(ind));
%     fitpart.Y(ind) = polyval(pfit.y,fitpart.t(ind));
%     fitpart.Z(ind) = polyval(pfit.z,fitpart.t(ind));
    
    % smoothing 
    fitpart.X(ind) = smooth(t(ind),part.X(ind),0.01,'rloess');
    fitpart.Y(ind) = smooth(t(ind),part.Y(ind),0.01,'rloess');
    fitpart.Z(ind) = smooth(t(ind),part.Z(ind),0.01,'rloess');

    % calculate velocity & accelearation from fitted coordinates
    fitpart.Vx(ind) = gradient(fitpart.X(ind),1/fps)./1000;
    fitpart.Vy(ind) = gradient(fitpart.Y(ind),1/fps)./1000;
    fitpart.Vz(ind) = gradient(fitpart.Z(ind),1/fps)./1000; % convert mm/s to m/s

    fitpart.Ax(ind) = gradient(fitpart.Vx(ind),1/fps);
    fitpart.Ay(ind) = gradient(fitpart.Vy(ind),1/fps);
    fitpart.Az(ind) = gradient(fitpart.Vz(ind),1/fps); % m/s^2

end

%% plot fit Vf
Vfinterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% vx
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Vfinterp.x = [Vfinterp.x; interp1(t(ind),fitpart.Vx(ind),tfinterp)];
    plot(t(ind),fitpart.Vx(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Vfinterp.x,"omitnan") + std(Vfinterp.x,"omitnan");
lowbound = mean(Vfinterp.x,"omitnan") - std(Vfinterp.x,"omitnan");
hm = plot(tfinterp,mean(Vfinterp.x,"omitnan"),'k--',LineWidth=2);

hs = patch([tfinterp fliplr(tfinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-0.35 0.35])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
title('$Velocity-MatlabSmoothFit$','interpreter','latex',FontWeight='bold',FontSize=18)
legend([hm,hs],'$\langle V_i \rangle$','$\langle V_i \rangle \pm \sigma (V_i)$','interpreter','latex',Location='south',FontSize=12);
ylabel('$V_x (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% vy
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Vfinterp.y = [Vfinterp.y; interp1(t(ind),fitpart.Vy(ind),tfinterp)];
    plot(t(ind),fitpart.Vy(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Vfinterp.y,"omitnan") + std(Vfinterp.y,"omitnan");
lowbound = mean(Vfinterp.y,"omitnan") - std(Vfinterp.y,"omitnan");
hm = plot(tfinterp,mean(Vfinterp.y,"omitnan"),'k--',LineWidth=2);

hs = patch([tfinterp fliplr(tfinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-0.7 0])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$V_y (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% vz
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Vfinterp.z = [Vfinterp.z; interp1(t(ind),fitpart.Vz(ind),tfinterp)];
    plot(t(ind),fitpart.Vz(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Vfinterp.z,"omitnan") + std(Vfinterp.z,"omitnan");
lowbound = mean(Vfinterp.z,"omitnan") - std(Vfinterp.z,"omitnan");
hm = plot(tfinterp,mean(Vfinterp.z,"omitnan"),'k--',LineWidth=2);

hs = patch([tfinterp fliplr(tfinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-0.35 0.35])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$V_z (m \cdot s^{-1})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim([-inf,0.35])

%
linkaxes(tf.Children,'x')
colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.')
set(col,'position',[0.42 0.34 0.2 0.02])

savefig('./Figures/Vf')
saveas(gcf,'./Figures/Vf','png')
saveas(gcf,'./Figures/Vf','pdf')

%% plot fit Af
Afinterp = struct('x',[],'y',[],'z',[]);
figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% vx
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Afinterp.x = [Afinterp.x; interp1(t(ind),fitpart.Ax(ind),tfinterp)];
    plot(t(ind),fitpart.Ax(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Afinterp.x,"omitnan") + std(Afinterp.x,"omitnan");
lowbound = mean(Afinterp.x,"omitnan") - std(Afinterp.x,"omitnan");
hm = plot(tfinterp,mean(Afinterp.x,"omitnan"),'k--',LineWidth=2);

hs = patch([tfinterp fliplr(tfinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-250 250])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
title('$Acce.-MatlabSmoothFit$','interpreter','latex',FontWeight='bold',FontSize=18)
legend([hm,hs],'$\langle A_i \rangle$','$\langle A_i \rangle \pm \sigma (A_i)$','interpreter','latex',Location='south',FontSize=12);
ylabel('$A_x (m \cdot s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% vy
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Afinterp.y = [Afinterp.y; interp1(t(ind),fitpart.Ay(ind),tfinterp)];
    plot(t(ind),fitpart.Ay(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Afinterp.y,"omitnan") + std(Afinterp.y,"omitnan");
lowbound = mean(Afinterp.y,"omitnan") - std(Afinterp.y,"omitnan");
hm = plot(tfinterp,mean(Afinterp.y,"omitnan"),'k--',LineWidth=2);

hs = patch([tfinterp fliplr(tfinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-250 250])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$A_y (m \cdot s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xticklabels([])
grid on
xlim([-inf,0.35])

% vz
nexttile;
for i = 1:numel(uni_exp)
    ind = find(nexp == uni_exp(i));
    Afinterp.z = [Afinterp.z; interp1(t(ind),fitpart.Az(ind),tfinterp)];
    plot(t(ind),fitpart.Az(ind),'-',Color=mycolormap(mycolorind(uni_exp(i)),:),LineWidth=2);hold on
end

upbound  = mean(Afinterp.z,"omitnan") + std(Afinterp.z,"omitnan");
lowbound = mean(Afinterp.z,"omitnan") - std(Afinterp.z,"omitnan");
hm = plot(tfinterp,mean(Afinterp.z,"omitnan"),'k--',LineWidth=2);

hs = patch([tfinterp fliplr(tfinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
ylim([-250 250])
patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
alpha(patchalpha1)
set(gca,FontSize=15)
ylabel('$A_z (m \cdot s^{-2})$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t(s)$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim([-inf,0.35])

%
linkaxes(tf.Children,'x')
colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
ylabel(col,'Num. Exp.')
set(col,'position',[0.42 0.34 0.2 0.02])

savefig('./Figures/Af')
saveas(gcf,'./Figures/Af','png')
saveas(gcf,'./Figures/Af','pdf')