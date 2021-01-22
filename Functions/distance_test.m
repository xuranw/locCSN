function [pval] = distance_test(D_mat, change, num_perm)
stat_rec = [];
n = size(D_mat, 1);

Z = zeros(n, 2);
Z(1:change, 1) = 1;
Z(change+1:end, 2) = -1;

counts = sum(abs(Z), 1);
count_mat = counts'*counts;

for i = 1:num_perm
    if i == 1
        ind = 1:n;
    else
        ind = randperm(n);
    end
    Dp = D_mat(ind, ind);
    test_stat = Z' * Dp*Z ./ (count_mat + 1e-6);
    stat_rec(i) = sum(test_stat(:));
end
pval = 1-mean(stat_rec(1) < stat_rec - 1e-6);
end