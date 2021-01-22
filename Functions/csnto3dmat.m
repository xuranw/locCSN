function [mat] = csnto3dmat(csn)
N = length(csn);
[N1, N2] = size(csn{1});
mat = zeros(N1, N2, N);
for i = 1:N
    mat(:, :, i) = csn{i};
end
end