function [overallSim, rxnSim] = FVAsimilarity(fva1, fva2)
v1mins=fva1(:,1);
v1maxs=fva1(:,2);
v2mins=fva2(:,1);
v2maxs=fva2(:,2);

rxnSim=[];
for counter=1:numel(v1mins)
    v1min=v1mins(counter)
    v1max=v1maxs(counter)
    v2min=v2mins(counter)
    v2max=v2maxs(counter)
    si=max(0,(min(v1max,v2max)-max(v1min,v2min)+eps)/(max(v1max,v2max)-min(v1min,v2min)+eps));
    rxnSim=[rxnSim; si];
end
overallSim = mean(rxnSim);
end