function csn = csn_origin_rec(chisq_mtx, data1, data2, wd_q, dev, md, iter, fuzzy)
tic;
if nargin == 4
    dev = false;
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
[G1, N] = size(data1);
[G2, ~] = size(data2);
csn = cell(1, N);
[I, J] = find(chisq_mtx);
L = length(I);
disp([num2str(L) ' pairs need calculation'])
csn_mat = zeros(L, N);
selected_gene = cell(1, 2);
selected_gene{1} = sort(unique(I)); selected_gene{2} = sort(unique(J));
grid_gene = cell(1, 2);
grid_gene{1} = zeros(G1, N); grid_gene{2} = zeros(G2, N);
gene_first_zero = cell(1, 2);
gene_first_zero{1} = zeros(1, G1); gene_first_zero{2} = zeros(1, G2);
for s = selected_gene{1}'
    gene_temp = data1(s,:);
    range_temp = range(gene_temp(gene_temp~=0));
    gene_cut = min(gene_temp(gene_temp~=0)): range_temp/20 : max(gene_temp); 
     
    if sum(gene_temp == 0)> 1
        gene_cut = [0, gene_cut];
        gene_first_zero{1}(s) = 1;
    end
    grid_gene{1}(s,:) = discretize(gene_temp, gene_cut);
end
for s = selected_gene{2}'
    gene_temp = data2(s,:);
    range_temp = range(gene_temp(gene_temp~=0));
    gene_cut = min(gene_temp(gene_temp~=0)): range_temp/20 : max(gene_temp); 
     
    if sum(gene_temp == 0)> 1
        gene_cut = [0, gene_cut];
        gene_first_zero{2}(s) = 1;
    end
    grid_gene{2}(s,:) = discretize(gene_temp, gene_cut);
end
if dev
    if fuzzy
        parfor m = 1:L
            i = I(m); j = J(m);
            gene1 = data1(i,:); gene2 = data2(j,:);data = [gene1; gene2];
    	
            grid1 = grid_gene{1}(i,:); grid2 = grid_gene{2}(j,:);
            grid_mix = grid1*21+grid2; 
            u_grid = unique(grid_mix); n_grid = zeros(1, length(u_grid));
            for s = 1:length(u_grid)
                n_grid(s) = sum(grid_mix == u_grid(s));
            end
            u_grid_mix = [floor((u_grid-1)/21); rem(u_grid-1, 21)+1; u_grid; n_grid];
            if gene_first_zero{1}(i)==1
                u_grid_mix = u_grid_mix(:, u_grid_mix(1,:)~=1);
            end
            if gene_first_zero{2}(j) == 1
                u_grid_mix = u_grid_mix(:, u_grid_mix(2,:)~=1);
            end
            cell_id = zeros(1, size(u_grid_mix, 2));
            cell_id_full = cell(1, size(u_grid_mix, 2));
            for t = 1:size(u_grid_mix, 2)
                cell_id_full{t} = find(grid_mix == u_grid_mix(3, t));
                cell_id(t) = randsample(cell_id_full{t}, 1);
            end
            [upper, lower] = upperlower_dev(gene1, gene2, wd_q, md, iter, cell_id);
            csn_temp = zeros(1, N);
            for t = 1:length(cell_id)
                k = cell_id(t);
                B = zeros(2, N);
                for l = 1:N
                    B(:, l) = data(:, l)<=upper(:, k) & data(:,l)>lower(:,k)&data(:,k);
                end
                B = sparse(B);
                a = sum(B, 2);
                temp = (B*B'*N-a*a')./sqrt((a*a').*((N-a)*(N-a)')/(N-1)+eps);
                csn_temp(:, cell_id_full{t}) = temp(1, 2);
            end
            csn_mat(m, :) = csn_temp;
            if rem(m, 10000) == 0
                disp([num2str(m) ' is done.'])
            end
        end
    else
        parfor m = 1:L
            i = I(m); j = J(m);
            gene1 = data1(i,:); gene2 = data2(j,:);data = [gene1; gene2];
            [upper, lower] = upperlower_dev(gene1, gene2, wd_q, md, iter);
            for k = 1:N
            if gene1(k)*gene2(k) > 0
                B = zeros(2, N);
                for l = 1:N
                    B(:, l) = data(:, l)<=upper(:,k) & data(:,l)>=lower(:,k)&data(:,k);
                end
                B = sparse(B);
                a = sum(B,2);
                temp = (B*B'*N-a*a')./sqrt((a*a').*((N-a)*(N-a)')/(N-1)+eps);
                csn_mat(m,k) = temp(1,2); 
            end
            end
        if rem(m, 10000) == 0
            disp([num2str(m) ' is done.'])
        end
        end
    end
end
for k = 1:N
    csn{k} = sparse(I, J, csn_mat(:, k), G1, G2);
end
toc;
end
