function csn = csn_origin_rec_loc(data1, data2, knn_index, wd_q, dev, md, iter)
tic;
if nargin == 4
    dev = false;
end
if nargin == 5 && dev
    md = 1; iter = true;
end
if nargin == 6 && dev
    iter = true;
end
[G1, N] = size(data1);
[G2, ~] = size(data2);
[nk, nc] = size(knn_index);
if nc~=N
    error('dimesion of data and knn not match')
end

csn = cell(1, N);
for i=1:N
    csn{i} = sparse(G1, G2);
end
parfor k = 1:N
    %tic;
    index_temp = knn_index(:, k);
    data1_sub = data1(:, index_temp); data2_sub = data2(:, index_temp);
    g1_index = find(data1_sub(:, 1)>0); L1_temp = length(g1_index);
    g2_index = find(data2_sub(:, 1)>0); L2_temp = length(g2_index);
    r = 1; 
    I = zeros(1, L1_temp*L2_temp); J = I; S = I;
    for i = 1:L1_temp
        for j = 1:L2_temp
            gi = g1_index(i); gj = g2_index(j);
            gene1 = data1_sub(gi, :); gene2 = data2_sub(gj, :);
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
            I(r) = gi; J(r) = gj; S(r) = temp(1,2);
            r = r+1;
            %csn{k}(gi, gj) = temp(1, 2);
        end
    end
    %toc;
    csn{k} = sparse(I, J, S, G1, G2);
end
toc;
end
