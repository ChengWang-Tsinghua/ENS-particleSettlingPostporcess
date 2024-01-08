function slipVeloCylind = slipVeloCylinderical(particle_part,tracers_part,neighborIdxAll,startTime,endTime,k)

%% terminal state range
Nframemax = 1e6;

kexp = 1;
if nargin>5
    kexp = k;
end
tstart = mod(startTime(kexp),Nframemax);
tend = mod(endTime(kexp),Nframemax);

slipVeloCylind = struct('rho',[],'theta',[],'z',[],'Urel',[],'Urelnorm',[]);
%% prepare idx for searching
idx_front_back = sort(vertcat(neighborIdxAll(tstart:tend).idx));

for i = 1:numel(idx_front_back)
    idxp = particle_part.Tf == tracers_part.Tf(idx_front_back(i));
    idxt = idx_front_back(i);
    
    Xp = [particle_part.Xf(idxp), particle_part.Yf(idxp), particle_part.Zf(idxp)];
    Vp = [particle_part.Vx(idxp), particle_part.Vy(idxp), particle_part.Vz(idxp)];
    Xf = [tracers_part.Xf(idxt), tracers_part.Yf(idxt), tracers_part.Zf(idxt)];
    Uf = [tracers_part.Vx(idxt), tracers_part.Vy(idxt), tracers_part.Vz(idxt)];
    normVp = Vp/norm(Vp);

    RM = rotationMatrix(Vp, [1,0,0]);
    
    [slipVeloCylind(i,:).rho, slipVeloCylind(i,:).theta, slipVeloCylind(i,:).z] = cartToCylind(Xp, Xf, RM);
    slipVeloCylind(i,:).Urel = Uf.*normVp-Vp;
end
