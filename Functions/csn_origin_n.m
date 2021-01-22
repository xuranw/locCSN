function [csn] = csn_origin_n(data_full, wd_q, dev, md, iter)
if nargin == 2
    dev = false;
end
if nargin == 3 && dev
    md = 1; iter = true;
end
if nargin == 4 && dev
    iter = true;
end
[n1, n2] = size(data_full);   %n1 genes, n2 cells
csn = cell(1, n2);
for i=1:n2
    csn{i} = zeros(n1);
end
tic;
if dev 
    for i = 1:n1-1
        for j = i+1:n1
            gene1 = data_full(i,:); gene2 = data_full(j,:); data = [gene1; gene2];
            [upper, lower] = upperlower_dev(gene1, gene2, wd_q, md, iter);
            for k = 1:n2
                if gene1(k)*gene2(k)>0
                    B = zeros(2, n2);
                    for l = 1:n2
                        B(:,l) = data(:,l) <= upper(:,k) & data(:,l)>=lower(:,k)&data(:,k);
                    end
                    B = sparse(B);
                    a = sum(B,2);
                    temp = (B*B'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1)+eps);
                    csn{k}(i, j) = temp(1, 2);
                    csn{k}(j, i) = csn{k}(i,j);
                end
            end
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
        temp(1:n1+1:n1^2) = 0;
        csn{k} = temp;
        %disp(['Cell ' num2str(k) ' is completed']);
    end
end
toc;
end
