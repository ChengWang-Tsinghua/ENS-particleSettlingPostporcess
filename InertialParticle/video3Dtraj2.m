function video3Dtraj2(particle_part,tracer_part,origin,trajlen,trajincrp,trajincrt,Rmin,Rmax,Vmin,Vmax)

mycolormap = mycolor('#0000B2','#FFFFFF','#B10000');

x0 = origin(1);
z0 = origin(3);
%%
videoFile = VideoWriter('traj3D_exp.mp4', 'MPEG-4');
videoFile.FrameRate = 30; % Adjust the frame rate as needed
open(videoFile);

tf = particle_part.Tf;
nframe = numel(tf);

for idx1 = 1:nframe
    disp(num2str(idx1/nframe))
    
    %% trajectory video
    tstartp = max(min(tf),tf(idx1)-trajlen);
    tendp = tf(idx1);
    
    % particle's traj
    trajp = tstartp:trajincrp:tendp;

    for j = 1:numel(trajp)
        idxp = find(tf == trajp(j));
        if ~isempty(idxp)
            scatter3(particle_part.Zf(idxp)-z0,particle_part.Xf(idxp)-x0,particle_part.Yf(idxp),150,particle_part.Vy(idxp),'filled');hold on
            plot3(particle_part.Zf(idxp)-z0,particle_part.Xf(idxp)-x0,particle_part.Yf(idxp),'go',MarkerSize=10);hold on
        end
    end
    
    % tracer's traj
    trajt = tstartp:trajincrt:tendp;
    for k  = 1:numel(trajt)
        idxt = find(tracer_part.Tf==trajt(k));
        if ~isempty(idxt)
            scatter3(tracer_part.Zf(idxt)-z0,tracer_part.Xf(idxt)-x0,tracer_part.Yf(idxt),5,tracer_part.Vy(idxt),'filled',LineWidth=0.05);hold on
        end
    end
   
    axis equal
    axis([-25 25 -25 25 -40 40 ])
    xlabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=15)
    ylabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=15)
    zlabel('$y(mm)$','interpreter','latex',FontWeight='bold',FontSize=15)
    colormap(mycolormap)
    col = colorbar(FontSize=12,Location='eastoutside');
    ylabel(col,'$V_y(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
    caxis([Vmin Vmax])
    %% find neighbors of current frame
    idx2 = find(tracer_part.Tf==tf(idx1));
    neighborAll = neighborIdx2(particle_part, tracer_part, idx1, idx2, Rmin, Rmax);
    
    plot3(tracer_part.Zf(neighborAll.idxfront)-z0,tracer_part.Xf(neighborAll.idxfront)-x0,tracer_part.Yf(neighborAll.idxfront),'yo',MarkerSize=8);
    plot3(tracer_part.Zf(neighborAll.idxback)-z0,tracer_part.Xf(neighborAll.idxback)-x0,tracer_part.Yf(neighborAll.idxback),'ko',MarkerSize=8);
    clear neighborAll
    hold off
    %%
    pause(0.5)
%     view([i 30])
    set(gcf, 'Color', [1 1 1]);
    set(gcf, 'Renderer', 'opengl');
    frameImage = getframe(gcf);
    frameout = './videoFrames/';
    if ~exist(frameout)
        mkdir(frameout)
    end
    fname = [frameout 'frame_' num2str(idx1) '.fig'];
    saveas(gcf,fname)
    writeVideo(videoFile, frameImage);
end
close(videoFile);


