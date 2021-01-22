function [csn] = csn_origin_loc(data_full, knn_index, wd_q, dev, md, iter)
if nargin == 3
    dev = false;
end
if nargin == 4 && dev
    md = 1; iter = true;
end
if nargin == 5 && dev
    iter = true;
end
[n1, n2] = size(data_full);   %n1 genes, n2 cells
[nk, nc] = size(knn_index);
if nc~=n2
    error('dimesion of data and knn not match')
end
csn = cell(1, n2);
for i=1:n2
    csn{i} = sparse(n1, n1);
end
tic;
parfor k = 1:n2
    %tic;
    index_temp = knn_index(:, k);
    data_sub = data_full(:, index_temp);
    g_index = find(data_sub(:, 1)>0); L_temp = length(g_index);
    I = zeros(1, L_temp*(L_temp-1)); J = I; S = I;
    r = 1;
    for i = 1:L_temp-1
        for j = i+1:L_temp
            gi = g_index(i); gj = g_index(j);
            gene1 = data_sub(gi, :); gene2 = data_sub(gj, :);
            data = [gene1; gene2];
            if dev
                [upper, lower] = upperlower_dev(gene1, gene2, wd_q, md, iter, 1);
            else
                [upper, lower] = upperlower(data, wd_q);
            end
            B = zeros(2, nk);
            for m = 1:nk
                B(:, m) = data(:,m) <= upper(:,1) & data(:,m)>=lower(:,1);
            end
            B = sparse(B);
            a = sum(B, 2);
            temp = (B*B'*nk-a*a')./sqrt((a*a').*((nk-a)*(nk-a)')/(nk-1)+eps);
            I(r) = gi; J(r) = gj; S(r) = temp(1, 2);
            r = r+1;
            I(r) = gj; J(r) = gi; S(r) = temp(1, 2);
            r = r+1;
            %csn{k}(gi, gj) = temp(1, 2);
            %csn{k}(gj, gi) = csn{k}(gi, gj);
        end
    end
    %toc;
    csn{k} = sparse(I, J, S, n1, n1);
end
toc;
end
