clear
clc
close all

fin = 'D:\Chronos_Footage\1.1g\output';
fout = 'D:\Chronos_Footage\1.1g\output';

%%
f = waitbar(0,'Please wait...Reading');
mag = load([fin '\particle_PTV.mat']);
tracer = load([fin '\tracer_STB.mat']);

dist.nt5 = [];
dist.nt10 = [];
dist.dist5sum = [];
dist.dist10sum= [];
dist.dist5mean = [];
dist.dist10mean= [];

waitbar(0, f, 'Please wait...');
for i = 1:size(mag.part.x,1)
    waitbar(i/size(mag.part.x,1), f, 'Please wait...');
    dist.dist5sum(i,1) = 0;
    dist.dist10sum(i,1) = 0;
    dist.nt5(i,1) = 0;
    dist.nt10(i,1) = 0;
    nframe = mag.part.nframe(i);
    Nexp = mag.part.Nexp(i);
    ind1 = find(tracer.part.nframe == nframe); % at the same frame
    ind2 = find(tracer.part.Nexp(ind1) == Nexp); % tracer number from the same experiment
    for j = 1:size(ind2)
        ind = ind2(j);
        d = sqrt((mag.part.x(i) - tracer.part.x(ind))^2 + (mag.part.y(i) - tracer.part.y(ind))^2 + (mag.part.z(i) - tracer.part.z(ind))^2);
        if d<5
            dist.nt5(i,1) = dist.nt5(i,1) + 1;
            dist.dist5sum(i,1) = dist.dist5sum(i,1) + d;
        elseif d<10
            dist.nt10(i,1) = dist.nt10(i,1) + 1;
            dist.dist10sum(i,1) = dist.dist10sum(i,1) + d;
        end
    end
    dist.dist5mean(i,1) = dist.dist5sum(i,1)/dist.nt5(i,1);
    dist.dist10mean(i,1) = dist.dist10sum(i,1)/dist.nt10(i,1);
end

for i = 1:max(mag.part.Nexp)
    ind3 = find(mag.part.Nexp == i);
    ind4 = find(dist.nt5(ind3) ~= 0);
    ind5 = find(dist.nt10(ind3) ~= 0);

    nt5(i,1) = mean(dist.nt5(ind3));
    nt10(i,1) = mean(dist.nt10(ind3));
    dist5(i,1) = mean(dist.dist5mean(ind4));
    dist10(i,1) = mean(dist.dist10mean(ind5));
end

waitbar(1, f, 'Done!');
pause(0.5)
close(f)


%%
load([fout '\NumInShell.mat'])

figure
subplot(2,1,1)
plot(dist.nt5)
title('Nt5')
subplot(2,1,2)
plot(dist.nt10)
title('Nt10')

figure
subplot(2,1,1)
plot(dist.dist5mean)
title('dist5mean')
subplot(2,1,2)
plot(dist.dist10mean)
title('dist10mean')

mean(nt5)
mean(nt10)

% save([fout '\NumInShell.mat'],'dist','dist5','dist10','nt5','nt10')

