clc;clear all
close all

% load

load('F:\MaltlabOutput\MultiCalib_1g_12V_0811_4kfps\particles\particle_PTV.mat')

fps = 4230;
nbins = 256;
dfit = 3;

%% fit 

Nexp = size(part,1);
Nframemax = 1e6;
NTrackmax = 1e9;
for nexp = 1:Nexp
    fit(nexp,:).T = part(nexp,:).T;
%     fit(nexp,:).T = (part(nexp,:).T-(nexp-1)*Nframemax)./fps;
    % fit x,y,z using polynomial
%     pfit(nexp,:).x= polyfit(part(nexp,:).t,part(nexp,:).x,dfit);
%     pfit(nexp,:).y= polyfit(part(nexp,:).t,part(nexp,:).y,dfit);
%     pfit(nexp,:).z= polyfit(part(nexp,:).t,part(nexp,:).z,dfit);
% 
%     fit(nexp,:).x = polyval(pfit(nexp,:).x,part(nexp,:).t);
%     fit(nexp,:).y = polyval(pfit(nexp,:).y,part(nexp,:).t);
%     fit(nexp,:).z = polyval(pfit(nexp,:).z,part(nexp,:).t);
    
    % smoothing 
    fit(nexp,:).X = smooth(fit(nexp,:).T,part(nexp,:).X,0.1,'rloess');
    fit(nexp,:).Y = smooth(fit(nexp,:).T,part(nexp,:).Y,0.1,'rloess');
    fit(nexp,:).Z = smooth(fit(nexp,:).T,part(nexp,:).Z,0.1,'rloess');

%     figure;
%     subplot(3,1,1);
%     plot(fit(nexp,:).T,part(nexp,:).X,'+');
%     axis tight;xlabel('t/s');ylabel('x/mm');hold on;
%     plot(fit(nexp,:).T,fit(nexp,:).X,LineWidth=2);
% 
%     subplot(3,1,2);
%     plot(part(nexp,:).T,part(nexp,:).Y,'+');
%     axis tight;xlabel('t/s');ylabel('y/mm');hold on;
%     plot(part(nexp,:).T,fit(nexp,:).Y,LineWidth=2);
% 
%     subplot(3,1,3);plot(fit(nexp,:).T,part(nexp,:).Z,'+');
%     axis tight;xlabel('t/s');ylabel('z/mm');hold on;
%     plot(fit(nexp,:).T,fit(nexp,:).Z,LineWidth=2);
%     
%     figure;
%     scatter3(part(nexp,:).X,part(nexp,:).Z,part(nexp,:).Y,'o')
%     hold on;
%     scatter3(fit(nexp,:).X,fit(nexp,:).Z,fit(nexp,:).Y)
%     axis equal

    % calculate velocity & accelearation from fitted coordinates
    fit(nexp,:).Vx = gradient(fit(nexp,:).X,1/fps)./1000;
    fit(nexp,:).Vy = gradient(fit(nexp,:).Y,1/fps)./1000;
    fit(nexp,:).Vz = gradient(fit(nexp,:).Z,1/fps)./1000; % convert mm/s to m/s

    fit(nexp,:).Ax = gradient(fit(nexp,:).Vx,1/fps);
    fit(nexp,:).Ay = gradient(fit(nexp,:).Vy,1/fps);
    fit(nexp,:).Az = gradient(fit(nexp,:).Vz,1/fps); % m/s^2

%     if max(abs(fit(nexp,:).Vx))>1 || max(abs(fit(nexp,:).Vy))>1 || max(abs(fit(nexp,:).Vz))>1
%         disp(num2str(nexp))
%     end

%     figure;
%     subplot(3,1,1);plot(fit(nexp,:).T,part(nexp,:).vx,'b+',fit(nexp,:).T,fit(nexp,:).Vx,'r-');
%     axis tight;xlabel('t/s');ylabel('Vx/[m/s]')
%     subplot(3,1,2);plot(fit(nexp,:).T,part(nexp,:).vy,'b+',fit(nexp,:).T,fit(nexp,:).Vy,'r-');
%     axis tight;xlabel('t/s');ylabel('Vy/[m/s]')
%     subplot(3,1,3);plot(fit(nexp,:).T,part(nexp,:).vz,'b+',fit(nexp,:).T,fit(nexp,:).Vz,'r-');
%     axis tight;xlabel('t/s');ylabel('Vz/[m/s]')
    
end

%% plot Velocity
figure
for nexp = 1:Nexp
%     plot(fit(nexp,:).T,part(nexp,:).vy,'+');
%     hold on
    plot(fit(nexp,:).T,fit(nexp,:).Vy,'-',LineWidth=2);hold on
end
xlim([0. 0.25])
ylim([-0.55 -0.1]);

%% stat
sumfit.Vx = [];
sumfit.Vy = [];
sumfit.Vz = [];
sumfit.Ax = [];
sumfit.Ay = [];
sumfit.Az = [];
for nexp = 1:Nexp
    sumfit.Vx = [sumfit.Vx;fit(nexp,:).Vx];
    sumfit.Vy = [sumfit.Vy;fit(nexp,:).Vy];
    sumfit.Vz = [sumfit.Vz;fit(nexp,:).Vz];
    sumfit.Ax = [sumfit.Ax;fit(nexp,:).Ax];
    sumfit.Ay = [sumfit.Ay;fit(nexp,:).Ay];
    sumfit.Az = [sumfit.Az;fit(nexp,:).Az];
end

%%
stats = stat_struct(sumfit);

gfitxV = linspace(-5,5,1028);
gfityV= normpdf(gfitxV,0,1);
gfitxA = linspace(-5,5,1028);
gfityA = normpdf(gfitxA,0,1);

%%
figure
[ipdf.Vxn,bin.Vxn] = FunPDF(stats.Vxn,nbins);
[ipdf.Vyn,bin.Vyn] = FunPDF(stats.Vyn,nbins);
[ipdf.Vzn,bin.Vzn] = FunPDF(stats.Vzn,nbins);
% [ipdf.Vn,bin.Vn] = FunPDF(stats.Vn,nbins);
semilogy(bin.Vxn,ipdf.Vxn,'ro',LineWidth=2);
hold on
semilogy(bin.Vyn,ipdf.Vyn,'go',LineWidth=2);
semilogy(bin.Vzn,ipdf.Vzn,'bo',LineWidth=2);
% semilogy(bin.Vn,ipdf.Vn,LineWidth=2);
semilogy(gfitxV,gfityV,'k',LineWidth=2)
legend('Vx','Vy-gravity','Vz','Gaussian')
xlabel('$(v_i-\langle v_i \rangle)/\sigma_{v_i}$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$p.d.f.$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
title('$Velocity$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlim([-5 5])
ylim([1e-3 1e0])


%% truncate
stats.Axn(abs(stats.Axn)>10)= [];
stats.Ayn(abs(stats.Ayn)>10)= [];
stats.Azn(abs(stats.Azn)>10)= [];
%%
figure
[ipdf.Axn,bin.Axn] = FunPDF(stats.Axn,nbins);
[ipdf.Ayn,bin.Ayn] = FunPDF(stats.Ayn,nbins);
[ipdf.Azn,bin.Azn] = FunPDF(stats.Azn,nbins);
% [ipdf.An,bin.An] = FunPDF(stats.An,nbins);
semilogy(bin.Axn,ipdf.Axn,'ro',LineWidth=2);
hold on
semilogy(bin.Ayn,ipdf.Ayn,'go',LineWidth=2);
semilogy(bin.Azn,ipdf.Azn,'bo',LineWidth=2);
% semilogy(bin.An,ipdf.An,LineWidth=2);
semilogy(gfitxA,gfityA,'k',LineWidth=1)
% xlabel('$(a_i-\langle a_i \rangle)/\sigma_{a_i} $',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlabel('$(a_i-\langle a_i \rangle)/\sigma_{a_i}$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
ylabel('$p.d.f.$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
legend('Ax','Ay-gravity','Az','Gaussian')
title('$Acceleration$',FontSize=20,FontName='Times New Roman',Interpreter='latex')
xlim([-6 6])
ylim([9e-4 2e1])