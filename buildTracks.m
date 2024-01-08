function [track,Ilong] = buildTracks(fpath,fname,ifParticle,frame0,long_thres)
%% Particle: load,remove fields,remove frames, build tracks, find long tracks
% load([fpathP_STB filesep 'STB_Particle']);
load([fpath fname]);

if ifParticle
    % remove fields
    fields = {'Vx','Vy','Vz','V','Ax','Ay','Az','A'};
    part = rmfield(part,fields);
    clear fields
end

% remove frames
% frame0 = 0;
part(1:frame0-1)=[];

% build track structure 
track = part2track(part);

% find long tracks
L = arrayfun(@(X)(numel(X.X)),track);
% long_thres = 10;
Ilong = find(L>=long_thres);

