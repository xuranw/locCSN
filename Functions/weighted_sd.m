function [sd_x] = weighted_sd(x, weight)
mu = sum(x.*weight)/sum(weight);
n = length(x);
sd_x = sqrt(sum((x - repmat(mu, 1, n)).^2.*weight)/sum(weight));
end