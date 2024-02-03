function dfolders = FunSubfolder(fin,vargin)
% return all subfolders

d =dir(fin);
% remove all files (isdir property is 0) and remove '.' and '..' 
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));

if nargin>1
    dfolders(cell2mat(arrayfun(@(X)(~endsWith(X.name,vargin)),dfolders,UniformOutput=false)))=[];
end
