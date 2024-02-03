function vec_pad = vertcat_pad(vec)

% if isrow(vec)
%     vec = vec';
% end

sizes = arrayfun(@(X)(numel(cell2mat(X))),vec);

vec_pad = [];
for i = 1:numel(vec)
    pad_nan = zeros(1,(max(sizes)-sizes(i)))*NaN;
    vec_pad = [vec_pad;[cell2mat(vec(i)),pad_nan]];
end