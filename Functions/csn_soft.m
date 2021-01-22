function [CSN, resample_cluster_store] = csn_soft(data_full, soft_c, wd_q, Nrep)
if nargin ==2
    wd_q = 0.1;
    Nrep = 100;
end
if nargin == 3
    Nrep = 100;
end

[N2, K] = size(soft_c);
CSN = cell(Nrep, N2);
resample_cluster_store = cell(1, Nrep);

for rep = 1:Nrep
tic;
resample_cluster = mnrnd(1, soft_c);
    resample_cluster_store{rep} = resample_cluster;
    csn = cell(1, N2);
    for cl = 1:K
        data = data_full(:, resample_cluster(:, cl)==1);
        [n1,n2] = size(data); %n1 genes, n2 cells.
        upper = zeros(n1,n2);
        lower = zeros(n1,n2);
        for i = 1 : n1
            [s1,s2] = sort(data(i,:));
            n3 = n2-sum(sign(s1));
            h = round(wd_q/2*sum(sign(s1)));
            k = 1;
            while k <= n2
                s = 0;
                while k+s+1 <= n2 && s1(k+s+1) == s1(k)
                    s = s+1;
                end
                if s >= h
                   upper(i,s2(k:k+s)) = data(i,s2(k));
                   lower(i,s2(k:k+s)) = data(i,s2(k));
                else
                   upper(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
                   lower(i,s2(k:k+s)) = data(i,s2(max(n3*(n3>h)+1,k-h)));
                end
                k = k+s+1;
            end
        end
        csn_temp = cell(1,n2);
        B = zeros(n1,n2);
        %p = -icdf('norm',alpha,0,1);
        for k = 1: n2  % for each cell
            for j = 1 : n2 % for each cell
                B(:,j) = data(:,j) <= upper(:,k) & data(:,j) >= lower(:,k);
            end
            a = sum(B,2);
            d = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
            d(1 : n1+1 : end) = 0;
            csn_temp{k} = d;
            %disp(['Cell ', num2str(k), ',cluster ', num2str(j), ' is completed']);
        end
        csn(resample_cluster(:, cl) == 1) = csn_temp;
    end
    CSN(rep, :) = csn;
    disp(['Rep', num2str(rep), ' is completed']);
toc;
end
end