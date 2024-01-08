function stats = stat_struct(part_struc)

% return stastatstistic of structure: mean,std,normalized

%% mean
stats.avgVx = mean(part_struc.Vx);
stats.avgVy = mean(part_struc.Vy);
stats.avgVz = mean(part_struc.Vz);
% stats.avgV = mean(part_struc.V);
stats.avgAx = mean(part_struc.Ax);
stats.avgAy = mean(part_struc.Ay);
stats.avgAz = mean(part_struc.Az);
% stats.avgA = mean(part_struc.A);

%% std
stats.stdVx = std(part_struc.Vx);
stats.stdVy = std(part_struc.Vy);
stats.stdVz = std(part_struc.Vz);
% stats.stdV = std(part_struc.V);
stats.stdAx = std(part_struc.Ax);
stats.stdAy = std(part_struc.Ay);
stats.stdAz = std(part_struc.Az);
% stats.stdA = std(part_struc.A);


%% original value

stats.Vx = part_struc;
stats.Vy = part_struc.Vy;
stats.Vz = part_struc.Vz;
% stats.V = part_struc.Vz;

stats.Ax = part_struc.Ax;
stats.Ay = part_struc.Ay;
stats.Az = part_struc.Az;
% stats.A = part_struc.Az;

%% normalized
stats.Vxn = (part_struc.Vx-stats.avgVx)/stats.stdVx;
stats.Vyn = (part_struc.Vy-stats.avgVy)/stats.stdVy;
stats.Vzn = (part_struc.Vz-stats.avgVz)/stats.stdVz;
% stats.Vn = (part_struc.Vz-stats.avgVz)/stats.stdV;

stats.Axn = (part_struc.Ax-stats.avgAx)/stats.stdAx;
stats.Ayn = (part_struc.Ay-stats.avgAy)/stats.stdAy;
stats.Azn = (part_struc.Az-stats.avgAz)/stats.stdAz;
% stats.An = (part_struc.Az-stats.avgAz)/stats.stdA;