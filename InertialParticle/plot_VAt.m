function [tinterp,tAverVyStart,tAverVyend,vy,vs,ay] = plot_VAt(tracklong_particle,startTime,endTime,axisrange,Fs,ifsave,fout)

%%
Nframemax = 1e6;

Nexp = max(size(tracklong_particle));


t = arrayfun(@(X)(mod(X.Tf,Nframemax)/Fs),tracklong_particle,'UniformOutput',false)';



%%
% tmin = min(t);
% tmax = max(t);

tmin = min(cell2mat(t));
tmax = max(cell2mat(t));
tinterp = linspace(tmin,tmax,1000);

tAverVyStart = tmin+startTime/Fs;
tAverVyend = tmin+endTime/Fs;

%%
mycolormap = mycolor('#063970','#eeeee4','#e28743');%('#063970','#eeeee4','#e28743')
mycolorind = floor(linspace(1,size(mycolormap,1),Nexp));
% patchcolor1 = mycolormap(size(mycolormap,1)/2,:);
patchcolor2 = mycolor('#a155b9');
% patchalpha1 = 0.4;
patchalpha2 = 0.3;

%% Velocity
Vinterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% vx
nexttile;
for kexp = 1:Nexp
    tk = t{kexp};
    Vx = tracklong_particle(kexp).Vx;
    Vinterp.x = [Vinterp.x; interp1(tk,Vx,tinterp)];
    plot(tk,Vx,'-',Color=mycolormap(mycolorind(kexp),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.x,"omitnan") + std(Vinterp.x,"omitnan");
lowbound = mean(Vinterp.x,"omitnan") - std(Vinterp.x,"omitnan");
hm = plot(tinterp,mean(Vinterp.x,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-0.15 0.15])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15)
legend([hm,hs],'$\langle V_i \rangle$','$\langle V_i \rangle \pm \sigma (V_i)$','interpreter','latex',Location='eastoutside',FontSize=12);
% title('$Velocity$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$V_x [mm/s]$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

% vy
nexttile;
for kexp = 1:Nexp
    tk = t{kexp};
    Vy = tracklong_particle(kexp).Vy;
    Vinterp.y = [Vinterp.y; interp1(tk,Vy,tinterp)];
    plot(tk,Vy,'-',Color=mycolormap(mycolorind(kexp),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.y,"omitnan") + std(Vinterp.y,"omitnan");
lowbound = mean(Vinterp.y,"omitnan") - std(Vinterp.y,"omitnan");
hm = plot(tinterp,mean(Vinterp.y,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-0.7 0])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15)
ylabel('$V_y [mm/s]$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

% vz
nexttile;
for kexp = 1:Nexp
    tk = t{kexp};
    Vz = tracklong_particle(kexp).Vz;
    Vinterp.z = [Vinterp.z; interp1(tk,Vz,tinterp)];
    plot(tk,Vz,'-',Color=mycolormap(mycolorind(kexp),:),LineWidth=2);hold on
end

upbound  = mean(Vinterp.z,"omitnan") + std(Vinterp.z,"omitnan");
lowbound = mean(Vinterp.z,"omitnan") - std(Vinterp.z,"omitnan");
hm = plot(tinterp,mean(Vinterp.z,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-0.15 0.15])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15)
ylabel('$V_z [mm/s]$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='eastoutside');
col.TickLabelInterpreter = "latex";
ylabel(col,'Num. Exp.')
% set(col,'position',[0.42 0.34 0.2 0.02])

%
% linkaxes(tf.Children,'x')

if ifsave ==1
    savefig_custom(fout,'V',8,6)
end
%% Acceleration
Ainterp = struct('x',[],'y',[],'z',[]);

figure
tf=tiledlayout(3,1,'TileSpacing','compact');

% ax
nexttile;
for kexp = 1:Nexp
    tk = t{kexp};
    Ax = tracklong_particle(kexp).Ax;
    Ainterp.x = [Ainterp.x; interp1(tk,Ax,tinterp)];
    plot(tk,Ax,'-',Color=mycolormap(mycolorind(kexp),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.x,"omitnan") + std(Ainterp.x,"omitnan");
lowbound = mean(Ainterp.x,"omitnan") - std(Ainterp.x,"omitnan");
hm = plot(tinterp,mean(Ainterp.x,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-250 250])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15)
legend([hm,hs],'$\langle A_i \rangle$','$\langle A_i \rangle \pm \sigma (A_i)$','interpreter','latex',Location='eastoutside',FontSize=12);
% title('$Acce.$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$A_x [mm/s^{2}]$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

% ay
nexttile;
for kexp = 1:Nexp
    tk = t{kexp};
    Ay = tracklong_particle(kexp).Ay;
    Ainterp.y = [Ainterp.y; interp1(tk,Ay,tinterp)];
    plot(tk,Ay,'-',Color=mycolormap(mycolorind(kexp),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.y,"omitnan") + std(Ainterp.y,"omitnan");
lowbound = mean(Ainterp.y,"omitnan") - std(Ainterp.y,"omitnan");
hm = plot(tinterp,mean(Ainterp.y,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-250 250])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15)
ylabel('$A_y [mm/s^{2}]$','interpreter','latex',FontWeight='bold',FontSize=18)
% xticklabels([])
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

% az
nexttile;
for kexp = 1:Nexp
    tk = t{kexp};
    Az = tracklong_particle(kexp).Az;
    Ainterp.z = [Ainterp.z; interp1(tk,Az,tinterp)];
    plot(tk,Az,'-',Color=mycolormap(mycolorind(kexp),:),LineWidth=2);hold on
end

upbound  = mean(Ainterp.z,"omitnan") + std(Ainterp.z,"omitnan");
lowbound = mean(Ainterp.z,"omitnan") - std(Ainterp.z,"omitnan");
hm = plot(tinterp,mean(Ainterp.z,"omitnan"),'k.',LineWidth=1);

hs = patch([tinterp fliplr(tinterp)], [upbound fliplr(lowbound)],patchcolor2,Edgecolor=patchcolor2,LineWidth=1);
alpha(patchalpha2)
% ylim([-250 250])
% patch([min(xlim) min(xlim) tmin tmin], [min(ylim) max(ylim) max(ylim) min(ylim)], patchcolor1)
% alpha(patchalpha1)
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15)
ylabel('$A_z [mm/s^{2}]$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
xlim(axisrange)

colormap(mycolormap)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',round([0 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='eastoutside');
col.TickLabelInterpreter = "latex";
ylabel(col,'Num. Exp.');
% set(col,'position',[0.42 0.34 0.2 0.02])

%
% linkaxes(tf.Children,'x')

if ifsave ==1
    savefig_custom(fout,'A',8,6)
end

%% estimate taup
vy = mean(Vinterp.y,"omitnan");
ay = mean(Ainterp.y,"omitnan");
figure;
subplot(2,1,1)
plot(tinterp,vy,'k-');hold on
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15);grid
xlim(axisrange)
ylabel('$\langle V_y(t) \rangle [mm/s]$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)

v0 = vy(1);
vs = mean(vy(tinterp<tAverVyend & tinterp>tAverVyStart),"omitnan");
vn = (vy-v0)/(vs-v0);
subplot(2,1,2);
plot(tinterp,vn,'k-');hold on
plot([tmin tmin],ylim,'k--',LineWidth=1)
plot([tAverVyStart tAverVyStart],ylim,'r-.',LineWidth=1)
plot([tAverVyend tAverVyend],ylim,'r-.',LineWidth=1)
set(gca,FontSize=15);grid
xlim(axisrange)
ylabel('$\frac{\langle V_y (t)\rangle-V_0)}{(V_s-V_0)}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t [s]$','interpreter','latex',FontWeight='bold',FontSize=18)
if ifsave ==1
    savefig_custom(fout,'Vs',8,4)
end
