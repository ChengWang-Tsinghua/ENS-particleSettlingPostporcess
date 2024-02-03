function plot_traj(tracklong_particle,origin,axisrange,ifsave,fout)

scale = 0.008;
Nexp = max(size(tracklong_particle));


x0 = origin(1);
z0 = origin(3);

mycolormap = mycolor('#B10000','#FFFFFF','#0000B2');
mycolormap2 = mycolor('#063970','#eeeee4','#e28743');
mycolorind = 1:floor(size(mycolormap,1)/Nexp):size(mycolormap,1);
%% all trajectories 3D plot
figure
for kexp = 1:Nexp
    scatter3(tracklong_particle(kexp).Zf-z0,tracklong_particle(kexp).Xf-x0,tracklong_particle(kexp).Yf,abs(tracklong_particle(kexp).Ay)*scale,tracklong_particle(kexp).Vy,'filled',LineWidth=2);hold on
end
axis equal;box;grid on;axis padded
set(gca,FontSize=15)
legend('$size \propto Acce. $','interpreter','latex',FontWeight='bold',FontSize=15)
% title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$y [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$x [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$z [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
legend(Location='northeast')
colormap(mycolormap);
col =colorbar;
col.TickLabelInterpreter = "latex";
ylabel(col,'$Vy(m/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
axis(axisrange)

if ifsave == 1
    savefig_custom(fout,'traj3D',8,7)
end
%% all trajectories 2D plot - xz projection 
figure
for kexp = 1:Nexp
    scatter(tracklong_particle(kexp).Xf-x0,tracklong_particle(kexp).Zf-z0,abs(tracklong_particle(kexp).Ay)*scale,tracklong_particle(kexp).Vy,'filled',LineWidth=2);hold on
end
scatter(0,0,'go',MarkerFaceColor='g',LineWidth=2)
axis equal;box;grid;axis padded
set(gca,FontSize=15)
legend('$size \propto Acce. $','interpreter','latex',FontWeight='bold',FontSize=15)
% title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap);
col =colorbar;
col.TickLabelInterpreter = "latex";
ylabel(col,'$Vy(m/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
axis(axisrange(1:4))

if ifsave == 1
    savefig_custom(fout,'traj2Dxz',8,7)
end

%% xz projection
figure
for kexp = 1:Nexp
    plot(tracklong_particle(kexp).Xf-x0,tracklong_particle(kexp).Zf-z0,Color=mycolormap2(mycolorind(kexp),:),LineWidth=2);hold on
end
scatter(0,0,'go',MarkerFaceColor='g',LineWidth=2)
axis equal;box on;grid on;axis padded
set(gca,FontSize=15)
% title(titlestring,'interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$z [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$x [mm]$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2)
col = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',ceil([0.01 0.25 0.5 0.75 1]*Nexp),FontSize=12,Location='north');
col.TickLabelInterpreter = "latex";
ylabel(col,'Num. Exp.')
set(col,'position',[0.65 0.2 0.2 0.02])
axis(axisrange(1:4))

if ifsave == 1
    savefig_custom(fout,'traj2Dxz2',8,7)
end