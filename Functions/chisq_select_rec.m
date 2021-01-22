function chisq_mtx = chisq_select_rec(c_data1, c_data2, alpha)
[G1, ~] = size(c_data1);
[G2, ~] = size(c_data2);
chisq_mtx = zeros(G1, G2);
tic;
parfor i = 1:G1
    j_temp = zeros(1, G2)
    for j = 1:G2
        [~, ~, p] = crosstab(c_data1(i, :), c_data2(j, :));
        if p < alpha
            j_temp(j) = 1;
        end
    end
    chisq_mtx(i, :) = j_temp;
    if rem(i, 100) == 0
        disp(['i = ' num2str(i) ' is completed']);
    end
end
toc;
[i, j, s] = find(chisq_mtx);
chisq_mtx = sparse(i, j, s, G1, G2);
end