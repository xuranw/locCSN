function [csn] = csn_origin(chisq_mtx, data_full, wd_q, dev, md, iter, fuzzy)
if nargin == 3
    dev = false;
end
if nargin == 4 && dev
    md = 1; iter = true; fuzzy = false;
end
if nargin == 5 && dev
    iter = true; fuzzy = false;
end
if nargin == 6 && dev
    fuzzy = false;
end
[n1, n2] = size(data_full);
csn = cell(1, n2);

tic;
[I, J] = find(chisq_mtx);
L = length(I);
disp([num2str(L) ' pairs need calculation'])
csn_mat = zeros(L, n2);
selected_gene = sort(unique([I; J]));
gene_mid = cell(1, n1);
grid_gene = zeros(n1, n2);
gene_first_zero = zeros(1, n1);
for s = selected_gene'
    gene_temp = data_full(s,:);
    range_temp = range(gene_temp(gene_temp~=0));
    gene_cut = min(gene_temp(gene_temp~=0)): range_temp/20 : max(gene_temp); 
    gene_mid{s} = gene_cut(1:20)+range_temp/40; 
    if sum(gene_temp == 0)> 1
        gene_cut = [0, gene_cut];
        gene_first_zero(s) = 1;
    end
    grid_gene(s,:) = discretize(gene_temp, gene_cut);
end
if dev
    if fuzzy
    parfor m = 1:L
        i = I(m); j = J(m);
        gene1 = data_full(i, :); gene2 = data_full(j,:); data = [gene1; gene2];
        grid1 = grid_gene(i,:); grid2 = grid_gene(j,:);
        grid_mix = grid1*21+grid2; 
        u_grid = unique(grid_mix); n_grid = zeros(1, length(u_grid));
        for s = 1:length(u_grid)
            n_grid(s) = sum(grid_mix == u_grid(s));
        end
        u_grid_mix = [floor((u_grid-1)/21); rem(u_grid-1, 21)+1; u_grid; n_grid];
        if gene_first_zero(i)==1
            u_grid_mix = u_grid_mix(:, u_grid_mix(1,:)~=1);
        end
        if gene_first_zero(j) == 1
            u_grid_mix = u_grid_mix(:, u_grid_mix(2,:)~=1);
        end
        cell_id = zeros(1, size(u_grid_mix, 2));
        cell_id_full = cell(1, size(u_grid_mix, 2));
        for t = 1:size(u_grid_mix, 2)
            cell_id_full{t} = find(grid_mix == u_grid_mix(3, t));
            cell_id(t) = randsample(cell_id_full{t}, 1);
        end
        [upper, lower] = upperlower_dev(gene1, gene2, wd_q, md, iter, cell_id);
        csn_temp = zeros(1, n2);
        for t = 1:length(cell_id)
            k = cell_id(t);
            B = zeros(2, n2);
            for l = 1:n2
                B(:, l) = data(:, l)<=upper(:, k) & data(:,l)>lower(:,k)&data(:,k);
            end
            B = sparse(B);
            a = sum(B, 2);
            temp = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
            csn_temp(cell_id_full{t}) = temp(1,2); 
        end
        csn_mat(m, :) = csn_temp;
    end
    for k = 1:n2
        csn{k} = sparse(I, J, csn_mat(:, k), n1, n1);
    end
    else
        parfor m = 1:L
        i = I(m); j = J(m);
        gene1 = data_full(i, :); gene2 = data_full(j,:); data = [gene1; gene2];
        [upper, lower] = upperlower_dev(gene1, gene2, wd_q, md, iter);
        for k = 1:n2
            if gene1(k)*gene2(k) > 0
                B = zeros(2, n2);
                for l = 1:n2
                    B(:, l) = data(:, l)<=upper(:,k) & data(:,l)>=lower(:,k)&data(:,k);
                end
                B = sparse(B);
                a = sum(B,2);
                temp = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
                csn_mat(m,k) = temp(1,2); 
            end
        end
        end 
    for k = 1:n2
        csn{k} = sparse(I, J, csn_mat(:, k), n1, n1);
    end
    end
else
    [upper, lower] = upperlower(data_full, wd_q);
    for k = 1:n2
        B = zeros(n1, n2);
        for j = 1 : n2
            B(:,j) = data_full(:,j) <= upper(:,k) & data_full(:,j) >= lower(:,k) & data_full(:,k);
        end
        B = sparse(B);
        a = sum(B,2);
        temp = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
        temp(1:n1+1:n1^2) = -100;
        csn{k} = temp;
        %disp(['Cell ' num2str(k) ' is completed']);
    end
end
toc;
end
