clc;clear all
close all

% load

% load('D:\Chronos_Footage\1.1g\output\tracer_STB.mat')
load('D:\Chronos_Footage0\1.1g\output\particle_PTV.mat')


%% plot tracks

MytrackID = 0;

NID = find(part.ID == MytrackID);
figure
scatter3(part.x(NID),part.z(NID),part.y(NID),[],part.ay(NID),LineWidth=2)
hold on
% plot3(part.x(NID),part.z(NID),part.y(NID),LineWidth=2)
h = colorbar;
ylabel(h,'$a_y[m/s^2] $',FontSize=15,FontName='Times New Roman',Interpreter='latex')
% caxis([-10 10])
axis equal
xlim([-25 25])
ylim([-25 25])
zlim([-40 40])
xlabel('$x/mm $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$z/mm$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
zlabel('$y(g)/mm $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$Tracks$',FontSize=20,FontName='Times New Roman',Interpreter='latex')

%% plot all tracks

maxTrackID = max(part.ID);
trackIDrange = maxTrackID/100;
figure
for n = 0:trackIDrange
    nid = find(part.ID == n);
    scatter3(part.x(nid),part.z(nid),part.y(nid),[],part.ay(nid),LineWidth=2)
    hold on
%     plot3(part.x(nid),part.z(nid),part.y(nid),LineWidth=2)
end
h = colorbar;
ylabel(h,'$a_y[m/s^2] $',FontSize=15,FontName='Times New Roman',Interpreter='latex')
% caxis([-10 10])
axis equal
xlim([-25 25])
ylim([-25 25])
zlim([-40 40])
xlabel('$x/mm $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$z/mm$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
zlabel('$y-gravity/mm $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$Tracks$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
