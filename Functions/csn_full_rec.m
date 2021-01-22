function csn_mtx = csn_full_rec(chisq_mtx, data1, data2, c_data1, c_data2, select_id, alpha, fuzzy, edges1, edges2)
tic;
if nargin < 8 || isempty(edges1) || isempty(edges2)
    fuzzy = false;
end
if nargin < 7 || isempty(alpha)
    alpha = 0.01;
end
[G1, N] = size(data1);
if nargin < 6 || isempty(select_id)
    select_id = 1:N;
end

[i, j]= find(chisq_mtx);
L = length(i);
disp([num2str(L) ' pairs need calculation'])

[G2, ~] = size(data2);
csn = zeros(L, N);
if fuzzy 
    parfor l = 1:L
        csn(l, :) = twogenecor(data1(i(l), :), data2(j(l), :), c_data1(i(l), :), c_data2(j(l), :), select_id, alpha, fuzzy, edges1{i(l)}, edges2{j(l)});
        if rem(l, 10000) == 0
            disp([num2str(l) ' is done.'])
        end
    end
else
    parfor l = 1:L
        csn(l, :) = twogenecor(data1(i(l), :), data2(j(l), :), c_data1(i(l), :), c_data2(j(l), :), select_id, alpha);
        if rem(l, 10000) == 0
            disp([num2str(l) ' is done.'])
        end
    end
end

K = length(select_id);
csn_mtx = cell(1, K);
for k = 1:K
    csn_mtx{k} = sparse(i, j, mean(csn(:, select_id(k)), 2)', G1, G2);
end
toc;
end