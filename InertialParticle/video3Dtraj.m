function video3Dtraj(particle_part,tracer_part,kexp,trajlen,trajincrp,trajincrt,Rmin,Rmax)



% trajlen = 100;
% trajincrp= 10;
% trajincrt = 1;
% kexp = 5;
% d = 1;
% [xs,ys,zs] = sphere;
% xs = xs*d/2;
% ys = ys*d/2;
% zs = zs*d/2;
% surf(xs,ys,zs);

mycolormap = mycolor('#0000B2','#FFFFFF','#B10000');

%%
Nframemax = 1e6;

%%
idxp = find(particle_part.Tf<kexp*Nframemax & particle_part.Tf>(kexp-1)*Nframemax);
idxt = find(tracer_part.Tf<kexp*Nframemax & tracer_part.Tf>(kexp-1)*Nframemax);

%%
ccodep = particle_part.Vy;
vpmax = max(ccodep(idxp));
vpmin = min(ccodep(idxp));

ccodet = tracer_part.Vy;
vtmax = max(ccodet(idxp));
vtmin = min(ccodet(idxp));

vmax = max(vpmax,vtmax);
vmin = min(vpmin,vtmin);

%%
videoFile = VideoWriter(['traj3D_exp' num2str(kexp) '.mp4'], 'MPEG-4');
videoFile.FrameRate = 30; % Adjust the frame rate as needed
open(videoFile);

for i = 1:numel(idxp)
    disp(num2str(i/numel(idxp)))
    
    idx1 = idxp(i);

    %% trajectory video
%     fig1 = figure;
    % particle's traj
    idx_tstartp = max(min(idxp),idx1-trajlen);
%     idx_tendp = min(max(idxp),idx1+trajlen); 
    idx_tendp =idx1;
    idx_trajp = idx_tstartp:trajincrp:idx_tendp;

    for j = 1:numel(idx_trajp)
%         cidx= max(1,ceil((ccodep(idx_trajp(j))-vpmin)/(vpmax-vpmin))*size(mycolormap,1));
%         colorp(:,:,1) = ones(size(xs,1))*mycolormap(cidx,1);
%         colorp(:,:,2) = ones(size(xs,1))*mycolormap(cidx,2); 
%         colorp(:,:,3) = ones(size(xs,1))*mycolormap(cidx,3); 
%         surf(part.Xf(idx_trajp(j))+xs,part.Zf(idx_trajp(j))+zs,part.Yf(idx_trajp(j))+ys,colorp);hold on
        scatter3(particle_part.Xf(idx_trajp(j)),particle_part.Zf(idx_trajp(j)),particle_part.Yf(idx_trajp(j)),40,particle_part.Vy((idx_trajp(j))),'filled');hold on
    end
    
    % tracer's traj
    tracer_start = particle_part.Tf(idx_tstartp);
    tracer_end = particle_part.Tf(idx_tendp);
    tracer_traj = tracer_start:trajincrt:tracer_end;
    for k  = 1:numel(tracer_traj)
        idx2 = find(tracer_part.Tf==tracer_traj(k));
        scatter3(tracer_part.Xf(idx2),tracer_part.Zf(idx2),tracer_part.Yf(idx2),5,tracer_part.Vy(idx2),'filled',LineWidth=0.05);hold on
    end
   
    axis equal
    axis([-25 25 -25 25 -40 40 ])
    xlabel('$x(mm)$','interpreter','latex',FontWeight='bold',FontSize=15)
    ylabel('$z(mm)$','interpreter','latex',FontWeight='bold',FontSize=15)
    zlabel('$y(mm)$','interpreter','latex',FontWeight='bold',FontSize=15)
    colormap(mycolormap)
    col = colorbar(FontSize=12,Location='eastoutside');
    ylabel(col,'$V_y(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
    caxis([vtmin vtmax])
    %% find neighbors of current frame
    idx3 = find(tracer_part.Tf==particle_part.Tf(idx1));
    neighborAll(i) = neighborInfos(particle_part,tracer_part,idx1,idx3,Rmin,Rmax);
    
    plot3(particle_part.Xf(idx1),particle_part.Zf(idx1),particle_part.Yf(idx1),'go');hold on
%     plot3(tracer.Xf(neighbor_global(i).idx),tracer.Zf(neighbor_global(i).idx),tracer.Yf(neighbor_global(i).idx),'ro');
    plot3(tracer_part.Xf(neighborAll(i).idxfront),tracer_part.Zf(neighborAll(i).idxfront),tracer_part.Yf(neighborAll(i).idxfront),'yo');
    plot3(tracer_part.Xf(neighborAll(i).idxback),tracer_part.Zf(neighborAll(i).idxback),tracer_part.Yf(neighborAll(i).idxback),'ko');
    
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
    fname = [frameout 'frame_' num2str(i) '.fig'];
    saveas(gcf,fname)
    writeVideo(videoFile, frameImage);
end
close(videoFile);


