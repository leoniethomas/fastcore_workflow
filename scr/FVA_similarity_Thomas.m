function res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs)
    res=[];
    for counter=1:numel(v1mins)
        v1min=v1mins(counter);
        v1max=v1maxs(counter);
        v2min=v2mins(counter);
        v2max=v2maxs(counter);
        si=max(0,(min(v1max,v2max)-max(v1min,v2min)+eps)/(max(v1max,v2max)-min(v1min,v2min)+eps));
        res=[res; si];
    end

    res = mean(res);
    
end