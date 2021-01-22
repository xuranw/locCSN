function [d] = deltacon(x1, x2, ee)
if nargin < 3 || isempty(ee)
   ee = 0.01;
end
N = size(x1, 1);
D1 = diag(sum(x1)); D2 = diag(sum(x2));

S1 = inv(eye(N) + ee^2*D1 - ee*x1);
S2 = inv(eye(N) + ee^2*D2 - ee*x2);

d = sqrt(sum(sum((sqrt(S1) - sqrt(S2)).^2)));
end
