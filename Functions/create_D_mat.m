function D_mat = create_D_mat(csn)
T = length(csn);
D_mat = zeros(T);
if issparse(csn{1})
    disp('sparse matrix')
    for i = 1:T-1
        tic;
        for j = i+1:T
            [idx, ~] = find(csn{i} - csn{j});
            [a, ~] = hist(idx, unique(idx));
            if isempty(a)
                D_mat(i, j) = 0; D_mat(j, i) = 0;
            else
            	D_mat(i, j) = max(a);
            	D_mat(j, i) = D_mat(i, j);
            end
        end
        toc;
    end
else
    disp('dense matrix')
    tic;
    for i = 1:T-1
        for j = i+1:T
            D_mat(i, j) = max(sum(abs(csn{i} - csn{j})));
            D_mat(j, i) = D_mat(i, j);
        end
    end
    toc;
end
end
