function [csn_comb] = csn_comb_cluster_new(csn, soft_c)
[n2, K] = size(soft_c);        %% K: clusters  n2: cells
scale = sqrt(soft_c(:,1).^2 + soft_c(:, 2).^2);
n1 = size(csn{1}, 1);
csn_comb = cell(1, n2);
for i = 1:n2
    csn_comb{i} = zeros(n1);
end
for k = 1:K
    for i = 1:n2
        csn_comb{i} = (csn_comb{i} + csn{i, k}*soft_c(i, k))/scale(i);
    end
end
end
