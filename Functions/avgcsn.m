function ctavg_mtx = avgcsn(chisq_mtx, data, c_data, cell_type_id, select_id, alpha, fuzzy, edges)
tic;
if nargin < 7 || isempty(edges)
    fuzzy = false;
end
if nargin < 6 || isempty(alpha)
    alpha = 0.01;
end
[G, N] = size(data);
if nargin < 5 || isempty(select_id)
    select_id = 1:N;
end

[i, j]= find(chisq_mtx);
L = length(i);
disp([num2str(L) ' pairs need calculation'])
csn = zeros(L, N);

if fuzzy
    parfor l = 1:L
        csn(l, :) = twogenecor(data(i(l), :), data(j(l), :), c_data(i(l), :), c_data(j(l), :), select_id, alpha, fuzzy, edges{i(l)}, edges{j(l)});
        if rem(l, 5000) == 0
            disp([num2str(l) ' is done.'])
        end
    end
else
    parfor l = 1:L
        csn(l, :) = twogenecor(data(i(l), :), data(j(l), :), c_data(i(l), :), c_data(j(l), :), select_id, alpha);
        if rem(l, 5000) == 0
            disp([num2str(l) ' is done.'])
        end
    end
end
select_ct_id = unique(cell_type_id(select_id));
K = length(select_ct_id);
ctavg_mtx = cell(1, K);
for k = 1:K
    ctavg_mtx{k} = sparse(i, j, mean(csn(:, cell_type_id == select_ct_id(k)), 2)', G, G);
end
toc;
end