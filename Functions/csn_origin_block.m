function [csn_mtx] = csn_origin_block(chisq_mtx, data, M, wd_q, dev, md, iter, fuzzy)
if nargin == 4
    dev = false; fuzzy = false;
end
if nargin == 5 && dev
    md = 1; iter = true; fuzzy = false;
end
if nargin == 6 && dev
    iter = true; fuzzy = false;
end
if nargin == 7 && dev
    fuzzy = false;
end

[G, K]= size(data);
n = ceil(G/M);
group_n = zeros(1, G);
for i = 1:n
   group_n((i-1)*M + 1: min(i*M, G)) = i;
end
csn_mtx = cell(1, K);
for k = 1:K
    csn_mtx{k} = sparse(G, G);
end

%tic;
for i = 1:n
    data_i = data(group_n == i, :);
    for j = i:n
        if i == j
            chisq_mtx_temp = chisq_mtx(group_n == i, group_n == i);
            csn_mtx_temp = csn_origin(chisq_mtx_temp, data_i, wd_q, dev, md, iter, fuzzy);
            parfor k = 1:K
                csn_mtx{k}(group_n == i, group_n == i) = csn_mtx_temp{k};
            end
        else
            data_j = data(group_n == j, :);
            chisq_mtx_temp = chisq_mtx(group_n == i, group_n == j);
            csn_mtx_temp = csn_origin_rec(chisq_mtx_temp, data_i, data_j, wd_q, dev, md, iter, fuzzy);
            parfor k = 1:K
                csn_mtx{k}(group_n == i, group_n == j) = csn_mtx_temp{k};
            end
        end
        disp(['block [' num2str(i) ',' num2str(j) '] finished!']);
    end
end
%toc;
end
