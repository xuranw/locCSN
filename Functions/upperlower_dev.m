function [upper, lower] = upperlower_dev(gene1, gene2, boxsize, md, iter, cell_id)
data = [gene1; gene2];
[n1,n2] = size(data);
if nargin < 6
    cell_id = 1:n2;
end
if nargin < 5
    iter = false; cell_id = 1:n2;
end
if nargin < 4
    iter = flase; md = 1; cell_id = 1:n2;
end
if size(cell_id, 1) > 1
    cell_id = cell_id';
end
[up_q, low_q] = upperlower(data, boxsize);
upper = zeros(n1,n2);
lower = zeros(n1,n2);
if iter
    maxiter = 10000; track_save = cell(1, 2);
    track_save{1} = zeros(n2, maxiter); track_save{2} = zeros(n2, maxiter);
    for k = cell_id
        if gene1(k)*gene2(k)>0
        d2_0 = md*sqrt(var(gene2(gene1<=up_q(1,k)&gene1>=low_q(1,k))));
        d1_0 = md*sqrt(var(gene1(gene2<=up_q(2,k)&gene2>=low_q(2,k))));
        d2_1 = md*sqrt(var(gene2(gene1<=gene1(k)+d1_0&gene1>=gene1(k)-d1_0)));
        d1_1 = md*sqrt(var(gene1(gene2<=gene2(k)+d2_0&gene2>=gene2(k)-d2_0)));
        count = 0;
        while sqrt((d2_0-d2_1)^2 + (d1_0-d1_1)^2) < 10^-5 && count <10000
            d2_0 = d2_1; d1_0 = d1_1;
            d2_1 = md*sqrt(var(gene2(gene1<=gene1(k)+d1_0&gene1>=gene1(k)-d1_0)));
            d1_1 = md*sqrt(var(gene1(gene2<=gene2(k)+d2_0&gene2>=gene2(k)-d2_0)));
            count = count + 1;
            track_save{1}(k, count) = d1_1; track_save{2}(k, count) = d2_1;
        end
        if count >= 10000
            save('./track_save.mat', 'track_save', '-v7.3');
            error('Error. \n Iteration at cell %d exceeds 10000', k);
        end
        upper(1, k) = gene1(k)+d1_1; upper(2, k) = gene2(k)+d2_1;
        lower(1, k) = gene1(k)-d1_1; lower(2, k) = gene2(k)-d2_1;
        end
    end
else
    for k = cell_id
        if gene1(k)*gene2(k)>0
        d2 = md*sqrt(var(gene2(gene1<=up_q(1,k)&gene1>=low_q(1,k))));
        d1 = md*sqrt(var(gene1(gene2<=up_q(2,k)&gene2>=low_q(2,k))));
        upper(1, k) = gene1(k)+d1; upper(2, k) = gene2(k)+d2;
        lower(1, k) = gene1(k)-d1; lower(2, k) = gene2(k)-d2;
        end
    end
end