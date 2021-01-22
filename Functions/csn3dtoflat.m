function [csn_flat] = csn3dtoflat(csn_3d)
[ng, ng1, nc] = size(csn_3d);
if ng~=ng1
    error('dimension of genes does not match!');
end
csn_flat = zeros(ng*(ng-1)/2, nc); k = 0;
for i = 1:ng-1
    for j = i+1:ng
        k = k+1;
        csn_flat(k,:) = csn_3d(i, j, :);
    end
end
end