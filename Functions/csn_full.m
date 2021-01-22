function csn_mtx = csn_full(chisq_mtx, data, c_data, select_id, alpha, fuzzy, edges)
tic;
if nargin < 6 || isempty(edges)
    fuzzy = false;
end
if nargin < 5 || isempty(alpha)
    alpha = 0.01;
end
[G, N] = size(data);
if nargin < 4 || isempty(select_id)
    select_id = 1:N;
end

[I, J]= find(chisq_mtx);
L = length(I);
disp([num2str(L) ' pairs need calculation'])
csn = zeros(L, N);

if fuzzy
    parfor l = 1:L
        csn(l, :) = twogenecor(data(I(l), :), data(J(l), :), c_data(I(l), :), c_data(J(l), :), select_id, alpha, fuzzy, edges{I(l)}, edges{J(l)});
        if rem(l, 10000) == 0
            disp([num2str(l) ' is done.'])
        end
    end
else
    parfor l = 1:L
        csn(l, :) = twogenecor(data(I(l), :), data(J(l), :), c_data(I(l), :), c_data(J(l), :), select_id, alpha);
        if rem(l, 10000) == 0
            disp([num2str(l) ' is done.'])
        end
    end
end
%select_ct_id = unique(cell_type_id(select_id));
K = length(select_id);
csn_mtx = cell(1, K);
for k = 1:K
    csn_mtx{k} = sparse(I, J, mean(csn(:, select_id(k)), 2)', G, G);
end
toc;
end