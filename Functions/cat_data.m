%% generate category data for test 
% data: gene * cell matrix
% soft clustering for each cell
% bin: which cells belong to which category 
function [c_data, edges] = cat_data(data, soft_c, bin)
if size(soft_c, 1)~=size(data, 2)
    error('size of cluster and data do not match.');
end
    

[G, N] = size(data);
K = size(soft_c, 2);
c_data = zeros(G, N);

index.zero = find(sum(data, 2) == 0); index.nonzero = find(sum(data, 2) > 0);
data_new = data(index.nonzero, :);
data_nonzero = data_new; 
if ~isempty(index.zero)
    c_data(index.zero, :) = ones(length(index.zero), N);
end

data_nonzero(data_nonzero == 0) = NaN;
data_minmax = [max(data, [], 2), min(data_nonzero, [], 2), min(data, [], 2)];
h = (data_minmax(:, 1) - data_minmax(:, 2))*bin;

index.one = index.nonzero(h == 0);
index.normal = index.nonzero(h > 0);

if ~isempty(index.one)
    for i = index.one'
        c_data(i, :) = ones(1, N);
        c_data(i, data(i, :) > 0) = 2;
    end
end

h = h(h>0); data_minmax = data_minmax(h>0, :);
data_nonzero = data_nonzero(h>0, :);
data_new = data_new(h>0,: );

nn = ceil(max(mean(data_nonzero, 2, 'omitnan')./h) + 0.5);
mm = ceil(max((data_minmax(:, 1) - mean(data_nonzero, 2, 'omitnan'))./h - 0.5));
edges = h*[-nn:0, 1:mm] + repmat(mean(data_nonzero, 2, 'omitnan') + h/2, 1, nn+mm+1);


for i = 1:length(h)
     di = discretize(data(index.normal(i), :), edges(i, :));
        c_data(index.normal(i), :) = di -min(di) + 1;
end
end