function [csn] = csn_soft_dev_rec(data1, data2, soft_c, upper1, lower1, upper2, lower2, wd_q, md, iter)
if nargin == 3
    wd_q = 0.1;
    [upper1, lower1] = upperlower_soft(data1, soft_c, wd_q);
    [upper2, lower2] = upperlower_soft(data2, soft_c, wd_q);
    md = 1; iter = false;
end
if nargin == 7
    wd_q = 0.1; md =1; iter = false;
end
if nargin == 8
    md = 1;iter = false;
end
if nargin == 9
    iter = false;
end
if sum(size(upper1))==0 || sum(size(lower1))==0
    [upper1, lower1] = upperlower_soft(data1, soft_c, wd_q);
end
if sum(size(upper2))==0 || sum(size(lower2))==0
    [upper2, lower2] = upperlower_soft(data2, soft_c, wd_q);
end
[~, K] = size(soft_c);        %% K: clusters
[n1, nc] = size(data1);       %% n1: genes in data1, nc: cells
[n2, ~] = size(data2);        %% n2: genes in data2
csn = cell(nc,K);
for cl = 1:K
    for k = 1:nc
        csn{k, cl} = zeros(n1, n2);
    end
end
for cl = 1:K
    n_cl = sum(soft_c(:,cl)); soft_cl = soft_c(:, cl);
    for i = 1:n1       % genes in data1
        tic;
        for j = 1:n2   % genes in data2
            nz_index = find(data1(i,:).*data2(j,:)>0 & soft_cl'>0);
            for k = nz_index  % k: index of cells, use non-zeros
                btw_i = data1(i, :) <= upper1{cl}(i, k) & data1(i, :) >= lower1{cl}(i, k);
                btw_j = data2(j, :) <= upper2{cl}(j, k) & data2(j, :) >= lower2{cl}(j, k);

                sdj_0 = md*weighted_sd(data2(j, btw_i), soft_cl(btw_i)');
                sdi_0 = md*weighted_sd(data1(i, btw_j), soft_cl(btw_j)');

                if iter
                    track = -ones(2, 10002); track(1, 1) = sdi_0; track(2, 1) = sdj_0;
                    btw_i = data1(i,:)<=data1(i,k)+sdi_0 & data1(i,:)>=data1(i,k)-sdi_0;
                    btw_j = data2(j,:)<=data2(j,k)+sdj_0 & data2(j,:)>=data2(j,k)-sdj_0;
                    sdj_1 = md*weighted_sd(data2(j, btw_i), soft_cl(btw_i)');
                    sdi_1 = md*weighted_sd(data1(i, btw_j), soft_cl(btw_j)');
                    track(1, 2) = sdi_1; track(2, 2) = sdj_1;
                    count = 0;
                    while sqrt((sdi_0-sdi_1)^2+(sdj_0-sdj_1)^2)<10^-8&&count<10000&&sdi_1>0&&sdj_1>0
                        sdj_0 = sdj_1; sdi_0 = sdi_1;
                        btw_i = data1(i,:)<=data1(i,k)+sdi_0 & data1(i,:)>=data1(i,k)-sdi_0;
                        btw_j = data2(j,:)<=data2(j,k)+sdj_0 & data2(j,:)>=data2(j,k)-sdj_0;
                        sdj_1 = md*weighted_sd(data2(j, btw_i), soft_cl(btw_i)');
                        sdi_1 = md*weighted_sd(data1(i, btw_j), soft_cl(btw_j)');
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
                nx = sum(soft_cl(data1(i, :)<data1(i, k)+sdi & data1(i, :)>data1(i, k)-sdi));
                ny = sum(soft_cl(data2(j, :)<data2(j, k)+sdj & data2(j, :)>data2(j, k)-sdj));
                nxy = sum(soft_cl(data1(i, :)<data1(i, k)+sdi & data1(i, :)>data1(i, k)-sdi & ...
                    data2(j, :)<data2(j, k)+sdj & data2(j, :)>data2(j, k)-sdj));
                rho_xy = nxy/n_cl - (nx/n_cl)*(ny/n_cl);
                sigma_xy = nx*ny*(n_cl-nx)*(n_cl-ny)/(n_cl^4*(n_cl-1));
                csn{k, cl}(i, j) = rho_xy/sqrt(sigma_xy);
            end
        end
        toc;
    end
    disp(['soft cluster ', num2str(cl)]);
end
end

