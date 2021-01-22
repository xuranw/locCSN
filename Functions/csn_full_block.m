function [csn_mtx, select_id] = csn_full_block(chisq_mtx, data, c_data, cell_type_id, M, select_type_id, alpha, fuzzy, edges)
if nargin < 8 || isempty(edges)
    fuzzy = false;
end
if nargin < 7 || isempty(alpha)
    alpha = 0.01;
end
if nargin < 6 || isempty(select_type_id)
    select_type_id = unique(cell_type_id);
end
select_id = find(ismember(cell_type_id, select_type_id));
[G, ~]= size(data);
n = ceil(G/M);
group_n = zeros(1, G);
for i = 1:n
   group_n((i-1)*M + 1: min(i*M, G)) = i; 
end
K = length(select_id);
csn_mtx = cell(1, K);
for k = 1:K
    csn_mtx{k} = sparse(G, G);
end
%tic;
for i = 1:n
    c_data_i = c_data(group_n == i, :);
    data_i = data(group_n == i, :);
    if fuzzy
        edges_i = edges(group_n == i);
    end
    for j = i:n
        if i == j
            chisq_mtx_temp = chisq_mtx((i-1)*M + 1: min(i*M, G), (i-1)*M + 1: min(i*M, G));
            if fuzzy
                csn_mtx_temp = csn_full(chisq_mtx_temp, data_i, c_data_i, select_id, alpha, fuzzy, edges_i);
            else
                csn_mtx_temp = csn_full(chisq_mtx_temp, data_i, c_data_i, select_id, alpha);
            end
            parfor k = 1:K
                csn_mtx{k}((i-1)*M + 1: min(i*M, G), (i-1)*M + 1: min(i*M, G)) = csn_mtx_temp{k};
            end
        else
            c_data_j = c_data(group_n == j, :);
            data_j = data(group_n == j, :);
            if fuzzy
                edges_j = edges(group_n == j);
            end
            chisq_mtx_temp = chisq_mtx((i-1)*M + 1: min(i*M, G), (j-1)*M + 1: min(j*M, G));
            if fuzzy
                csn_mtx_temp = csn_full_rec(chisq_mtx_temp, data_i, data_j, c_data_i, c_data_j, select_id, alpha, fuzzy, edges_i, edges_j);
            else
                csn_mtx_temp = csn_full_rec(chisq_mtx_temp, data_i, data_j, c_data_i, c_data_j, select_id, alpha);
            end
            parfor k = 1:K
                csn_mtx{k}((i-1)*M + 1: min(i*M, G), (j-1)*M + 1: min(j*M, G)) = csn_mtx_temp{k};
            end
        end
        disp(['block [' num2str(i) ',' num2str(j) '] finished!']);
    end
end
%toc;
end