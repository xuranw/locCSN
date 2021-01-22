%% Check the correlation between two genes for each soft cluster

function [corlist] = cor_check(data, soft_c)
if size(data, 2) ~= size(soft_c, 1)
    error('lengths of two genes differ! Please check.');
end
G = size(data, 1); 
K = size(soft_c, 2);

corlist = zeros((G-1)*G/2, K);

for k = 1:K
    counter = 1;
    tic;
    for g = 1:G-1
        for j = g+1:G
            xtemp = data([g, j], :);
            p = kde(xtemp, 'rot', soft_c(:, k)');
            %p1 = marginal(p, [1]);
            %p2 = marginal(p, [2]);
            hp = hist(p);
            hp = hp/sum(sum(hp));
            [N, M] = size(hp);
            hp1 = sum(hp, 1); hp2 = sum(hp, 2);
            E1 = sum((1:N).*hp1); E2 = sum((1:M).*hp2');
            V1 = sum((1:N).^2.*hp1) - E1^2; 
            V2 = sum((1:M).^2.*hp2') - E2^2;
            Cov12 = sum(sum((1:N)'*(1:M).*hp)) - E1*E2;
            corlist(counter, k) = Cov12/sqrt(V1*V2); 
            counter = counter + 1;
        end
    end
    toc;
end

end