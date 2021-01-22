% Function ctgrdata is used for catograzie data set.
%It generates c_data which is the catograised data. with bin size bin

function [c_data, edges] = ctgrdata(data, bin, mid_grid)
if nargin < 3
  mid_grid = false;
end
[G, N] = size(data);
c_data = zeros(G, N);

index.zero = find(sum(data, 2) == 0);    % index of gene that are all zero
index.nonzero = find(sum(data, 2) > 0);  % index of gene that expressed
data_new = data(index.nonzero, :);       % expression of gene that expressed
data_nonzero = data_new; 
if ~isempty(index.zero)
    c_data(index.zero, :) = ones(length(index.zero), N);
end

data_nonzero(data_nonzero == 0) = NaN;
data_minmax = [max(data_new, [], 2), min(data_nonzero, [], 2), min(data_new, [], 2)];
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

if mid_grid
    Edges = cell(1, G);
    for i = 1:length(h)
        di = discretize(data_new(i,:), edges(i,:));
        lower = max(0, edges(i, unique(di)));
        upper = min(data_minmax(i, 1), edges(i, unique(di)+1)); 
        c_data(index.normal(i), :) = di - min(di) + 1;
        Edges{index.normal(i)} = [unique(c_data(index.normal(i),:)); (lower + upper)/2];
    end
    edges = Edges;
else
    for i = 1:length(h)
        di = discretize(data(index.normal(i), :), edges(i, :));
        c_data(index.normal(i), :) = di -min(di) + 1;
    end
end
end