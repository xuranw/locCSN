function [res] = mvncsn(mu, sigma, Dx, Dy, N)
if nargin < 5
    N = 1000;
end
%y = mvnpdf(x, mu, sigma);
sigma1 = sigma(1, 1); sigma2 = sigma(2, 2);
mu1 = mu(1); mu2 = mu(2); 

x_grid = mu1-3.5*sqrt(sigma1):0.01*sqrt(sigma1):mu1+3.5*sqrt(sigma1);
y_grid = mu2-3.5*sqrt(sigma2):0.01*sqrt(sigma2):mu2+3.5*sqrt(sigma2);

p_grid = zeros(length(x_grid), length(y_grid));
for i = 1:length(x_grid)
    for j = 1:length(y_grid)
        p_grid(i,j) = mvnpdf([x_grid(i), y_grid(j)], mu, sigma);
    end
end
p_grid = p_grid/sum(sum(p_grid));
grid.x = mu1-3*sqrt(sigma1):0.3*sqrt(sigma1):mu1+3*sqrt(sigma1);
grid.y = mu2-3*sqrt(sigma2):0.3*sqrt(sigma2):mu2+3*sqrt(sigma2);
test_stat = cell(length(grid.x), length(grid.y));
for xi = 1:length(grid.x)
    tic;
    for yi = 1:length(grid.y)
        temp = zeros(length(Dx), length(Dy));
        x1 = grid.x(xi); x2 = grid.y(yi);
        for i = 1:length(Dx)
            for j = 1:length(Dy)
                dx = Dx(i); dy = Dy(j);
                xid = find(x_grid<x1+dx & x_grid>x1-dx);
                yid = find(y_grid<x2+dy & y_grid>x2-dy);
                px = sum(sum(p_grid(xid, :))); py = sum(sum(p_grid(:, yid)));
                pxy = sum(sum(p_grid(xid, yid)));
                rho_xy = pxy - px*py;
                sd_xy = sqrt(px*py*(1-px)*(1-py)/N);
                temp(i, j) = rho_xy/sd_xy;
            end
        end
        test_stat{xi, yi} = temp;
        disp(['grid ', num2str(xi), ' and ', num2str(yi), ' is done'])
    end
toc;
end
res.test_stat = test_stat;
res.grid = grid;
end
