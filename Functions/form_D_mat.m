function D_mat = create_D_mat(csn)
T = length(csn);
D_mat = zeros(T);
if T < 300
tic;
for i = 1:T-1
    for j = i+1:T
        D_mat(i, j) = max(sum(abs(csn{i} - csn{j})));
        D_mat(j, i) = D_mat(i, j);
    end
end
toc;
else
tic;
K = (T-1)*T/2; d_arr = zeros(1, K);
parfor k=1:K
	i = ceil((1+sqrt(1+8*k))/2); j = k - (i-2)*(i-1)/2;
	d_arr(k) = max(sum(abs(csn{i} - csn{j})));
end
for k = 1:K
	i = ceil((1+sqrt(1+8*k))/2); j = k - (i-2)*(i-1)/2;
	D_mat(i, j) = d_arr(k);
end
D_mat = D_mat + D_mat';
toc;
end
end
