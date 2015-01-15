function [gaps, Ws, Wref] = calc_gaps(data, g)
    [N, d] = size(data);
    
    ks = cellfun(@(gm)gm.NComponents, g);       % K
    
    % Gap statistic
    gaps = zeros(length(g), 1);
    
    assgns = cell(length(g), 1);
    Ws = zeros(length(g), 1);
    Wref = zeros(length(g), 1);
    for i = 1:length(g)
        assgns{i} = cluster(g{i}, data);
        
        W = 0;
        for k = 1:max(assgns{i})
            cdata = data(assgns{i} == k, :);
            W = W + sum(sum((cdata - repmat(mean(cdata), size(cdata, 1), 1)).^2, 2));
        end
        Ws(i) = log(W);
        
        Wref(i) = log(d*N)-(2/d)*log(ks(i));
        
        gaps(i) = Wref(i)-Ws(i);
    end
end