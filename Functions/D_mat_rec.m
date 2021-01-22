function D_mat = D_mat_rec(csn1, csn2)
T1 = length(csn1); T2 = length(csn2);
D_mat = zeros(T1, T2);
if max(T1, T2) < 200
    tic;
    for i = 1:T1
        for j = 1:T2
            D_mat(i, j) = max(sum(abs(csn1{i} - csn2{j})));
        end
    end
    toc;
else
    tic;
    K = T1*T2; d_arr = zeros(1, K);
    parfor k = 1:K
        i = floor((k-1)/T2) + 1; j = k - (i-1)*T2;
        d_arr(k) = max(sum(abs(csn1{i} - csn2{j})));
    end
    for k = 1:K
        i = floor((k-1)/T2) + 1; j = k - (i-1)*T2;
        D_mat(i, j) = d_arr(k);
    end
    toc;
end
end
