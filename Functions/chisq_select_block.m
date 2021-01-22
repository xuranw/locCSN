function chisq_mtx = chisq_select_block(c_data, M, alpha)
tic;
[G, ~] = size(c_data);
chisq_mtx = zeros(G, G);
n = ceil(G/M);
group_n = zeros(1, G);
for i = 1:n
   group_n((i-1)*M + 1: min(i*M, G)) = i; 
end
for i = 1:n
    c_data_i = c_data(group_n == i, :);
    %[m_b_i, ~] = size(c_data_i);
    for j = i:n
        if i == j
            chisq_mtx_temp = chisq_select(c_data_i, alpha);
            chisq_mtx((i-1)*M + 1: min(i*M, G), (i-1)*M + 1: min(i*M, G)) = chisq_mtx_temp;
        else
            c_data_j = c_data(group_n == j, :);
            chisq_mtx_temp = chisq_select_rec(c_data_i,c_data_j, alpha);
            chisq_mtx((i-1)*M + 1: min(i*M, G), (j-1)*M + 1: min(j*M, G)) = chisq_mtx_temp;
        end
        disp(['block [' num2str(i) ',' num2str(j) '] is done!']);
    end
end
toc;
end
