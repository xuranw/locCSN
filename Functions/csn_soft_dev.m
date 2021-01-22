function [csn] = csn_soft_dev(data_full, soft_c, upper, lower, wd_q, md, iter)
if nargin ==2
    wd_q = 0.1;
    [upper, lower] = upperlower_soft(data_full, soft_c, wd_q); 
    md = 1; iter = false;
end
if nargin == 4
    wd_q = 0.1; md =1; iter = false;
end
if nargin == 5
    md = 1;iter = false;
end
if nargin == 6
    iter = false;
end
if sum(size(upper))==0 || sum(size(lower))==0
    [upper, lower] = upperlower_soft(data_full, soft_c, wd_q);
end
[~, K] = size(soft_c);        %% K: clusters
[n1, n2] = size(data_full);   %% n1: genes, n2: cells
csn = cell(n2,K);
for cl = 1:K
    for k = 1:n2
        csn{k, cl} = zeros(n1);
    end
end
for cl = 1:K
    n_cl = sum(soft_c(:,cl)); soft_cl = soft_c(:, cl);
    for i = 1:n1-1
        tic;
        for j = i+1:n1
            nz_index = find(data_full(i,:).*data_full(j,:)>0 & soft_cl'>0);
            for k = nz_index  % k: index of cells
                btw_i = data_full(i, :) <= upper{cl}(i, k) & data_full(i, :) >= lower{cl}(i, k);
                btw_j = data_full(j, :) <= upper{cl}(j, k) & data_full(j, :) >= lower{cl}(j, k);
                
                sdj_0 = md*weighted_sd(data_full(j, btw_i), soft_cl(btw_i)');
                sdi_0 = md*weighted_sd(data_full(i, btw_j), soft_cl(btw_j)');
                
                if iter
                    track = -ones(2, 10002); track(1, 1) = sdi_0; track(2, 1) = sdj_0;
                    btw_i = data_full(i,:)<=data_full(i,k)+sdi_0 & data_full(i,:)>=data_full(i,k)-sdi_0;
                    btw_j = data_full(j,:)<=data_full(j,k)+sdj_0 & data_full(j,:)>=data_full(j,k)-sdj_0;
                    sdj_1 = md*weighted_sd(data_full(j, btw_i), soft_cl(btw_i)');
                    sdi_1 = md*weighted_sd(data_full(i, btw_j), soft_cl(btw_j)');
                    track(1, 2) = sdi_1; track(2, 2) = sdj_1;
                    count = 0;
                    while sqrt((sdi_0-sdi_1)^2+(sdj_0-sdj_1)^2)<10^-8&&count<10000&&sdi_1>0&&sdj_1>0
                        sdj_0 = sdj_1; sdi_0 = sdi_1;
                        btw_i = data_full(i,:)<=data_full(i,k)+sdi_0 & data_full(i,:)>=data_full(i,k)-sdi_0;
                        btw_j = data_full(j,:)<=data_full(j,k)+sdj_0 & data_full(j,:)>=data_full(j,k)-sdj_0;
                        sdj_1 = md*weighted_sd(data_full(j, btw_i), soft_cl(btw_i)');
                        sdi_1 = md*weighted_sd(data_full(i, btw_j), soft_cl(btw_j)');
                        count = count + 1;
                        track(1,count+2) = sdi_1; track(2, count+2) = sdj_1;
                    end
                    if count>=10000
                        error('Error. \nIteration of Cluster %d at gene %d and gene %d at cell %d have exceed 10000 ',cl, i, j, k)
                    end
                    sdj = sdj_1; sdi = sdi_1;
                else
                    sdj = sdj_0; sdi = sdi_0;
                end
                nx = sum(soft_cl(data_full(i, :)<data_full(i, k)+sdi & data_full(i, :)>data_full(i, k)-sdi));
                ny = sum(soft_cl(data_full(j, :)<data_full(j, k)+sdj & data_full(j, :)>data_full(j, k)-sdj));
                nxy = sum(soft_cl(data_full(i, :)<data_full(i, k)+sdi & data_full(i, :)>data_full(i, k)-sdi & ...
                    data_full(j, :)<data_full(j, k)+sdj & data_full(j, :)>data_full(j, k)-sdj));
                rho_xy = nxy/n_cl - (nx/n_cl)*(ny/n_cl);
                sigma_xy = nx*ny*(n_cl-nx)*(n_cl-ny)/(n_cl^4*(n_cl-1));
                csn{k, cl}(i, j) = rho_xy/sqrt(sigma_xy);
                csn{k, cl}(j, i) = csn{k, cl}(i, j);
            end
        end
    toc;
    end
    disp(['soft cluster ', num2str(cl)]);
end
end

