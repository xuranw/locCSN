function D_mat = D_mat(csn)
T = length(csn);
D_mat = zeros(T);
if T < 200
tic;
for i = 1:T-1
    for j = i+1:T
        D_mat(i, j) = max(sum(abs(csn{i} - csn{j})));
    end
end
toc;
else
tic;
K = (T-1)*T/2; d_arr = zeros(1, K);
parfor k=1:K
	j = ceil((1+sqrt(1+8*k))/2); i = k - (j-2)*(j-1)/2;
	d_arr(k) = max(sum(abs(csn{i} - csn{j})));
end
for k = 1:K
	j = ceil((1+sqrt(1+8*k))/2); i = k - (j-2)*(j-1)/2;
	D_mat(i, j) = d_arr(k);
end
toc;
end
end
