function chisq_mtx = chisq_select(c_data, alpha)
[G, ~] = size(c_data);
chisq_mtx = zeros(G, G);
tic;
parfor i = 1:G-1
    j_temp = zeros(1, G)
    for j = i+1:G
        [~, ~, p] = crosstab(c_data(i, :), c_data(j, :));
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
chisq_mtx = sparse(i, j, s, G, G);
end