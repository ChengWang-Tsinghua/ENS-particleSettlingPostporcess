function averES = averEulerianStats(EulerStats)

fields0 = fieldnames(EulerStats);
fields1 = {'filtL','filtW','r','S2x','S2y','S2z','Sau','Saulong','Ruur','Ruu','PSDk','PSD'};
fields2 ={'Splong','SplongAbs'};
fields3 = setdiff(fields0,[fields1,fields2]);

for i = 1:numel(fields1)
    averES.(fields1{i}) = mean(EulerStats.(fields1{i}),1,"omitnan");
end

for i = 1:numel(fields2)
    for j = 1:size(EulerStats.(fields2{i}),2)
        averES.(fields2{i}){:,j} = mean(cell2mat(arrayfun(@(X)(cell2mat(X.(fields2{i})(:,j)')'),EulerStats,UniformOutput=false)),1,"omitnan");
    end
end

for i = 1:numel(fields3)
    averES.(fields3{i}) = EulerStats.(fields3{i});
end