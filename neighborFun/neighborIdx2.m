function [neighborIdxAll] = neighborIdx2(particle_part,tracer_part,idx1,idx2,Rmin,Rmax)

particle.xf = particle_part.Xf(idx1);
particle.yf = particle_part.Yf(idx1);
particle.zf = particle_part.Zf(idx1);
particle.vx = particle_part.Vx(idx1);
particle.vy = particle_part.Vy(idx1);
particle.vz = particle_part.Vz(idx1);

tracer.xf = tracer_part.Xf(idx2);
tracer.yf = tracer_part.Yf(idx2);
tracer.zf = tracer_part.Zf(idx2);

%% neighboring tracers

neighborIdxAll = struct('idx',[],'d',[]);
neighborIdxAll.Rmin = Rmin;
neighborIdxAll.Rmax = Rmax;


for i = 1:numel(tracer.xf)
    d(i,:) = sqrt((particle.xf-tracer.xf(i))^2+(particle.yf-tracer.yf(i))^2+(particle.zf-tracer.zf(i))^2);
end

%% global neighboring
idx01 = find(d>Rmin & d<Rmax);
neighborIdxAll.idx = idx2(idx01);
neighborIdxAll.d = d(idx01);


%% front and back 
rpt.x = tracer.xf(idx01) - particle.xf;
rpt.y = tracer.yf(idx01) - particle.yf;
rpt.z = tracer.zf(idx01) - particle.zf;

rv = particle.vx*rpt.x + particle.vy*rpt.y + particle.vz*rpt.z;
neighborIdxAll.idxfront = neighborIdxAll.idx(rv>0);
neighborIdxAll.idxback = neighborIdxAll.idx(rv<0);





