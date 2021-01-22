function [D_mat] = create_D_mat_fro(csn)
T = length(csn);
D_mat = zeros(T);
tic;
for i = 1:T-1
    for j = i+1:T
        D_mat(i, j) = norm(csn{i}-csn{j}, 'fro');
        D_mat(j, i) = D_mat(i, j);
    end
end
toc;
end