function [upper, lower] = upperlower_soft(data_full, soft_c, wd_q)
% soft_c: cells * clusters
[n2, K] = size(soft_c);

F_c = soft_c./repmat(sum(soft_c,1), n2, 1);
upper = cell(1, K); lower = cell(1, K);
for cl = 1:K
    fc = F_c(:, cl);
    [n1, n2] = size(data_full);
    upper{cl} = zeros(n1, n2); lower{cl} = zeros(n1, n2);
    for i = 1:n1
        [s1, s2] = sort(data_full(i,:));
        n3 = sum(s1 ==0);
        k = 1;
        while k<= n2
            s = 0;
            while k+s+1<=n2 && s1(k+s+1)==s1(k)
                s = s+1;
            end
            if sum(fc(s2(k:k+s)))>=wd_q/2
                upper{cl}(i, s2(k:k+s)) = data_full(i, s2(k));
                lower{cl}(i, s2(k:k+s)) = data_full(i, s2(k));
            else
                h=1;
                while h+k+s<=n2 && sum(fc(s2(k:k+s+h)))<wd_q/2
                    h = h+1;
                end
                upper{cl}(i, s2(k:k+s)) = data_full(i, s2(min(n2, k+h+s)));
                h = 1;
                while k-h>0 && sum(fc(s2(k-h:k)))<wd_q/2
                    h = h+1;
                end
                lower{cl}(i, s2(k:k+s)) = data_full(i, s2(max(n3*(n3>h)+1, k-h)));
            end
            k = k+s+1;
        end
        %disp(['soft cluster ', num2str(cl), ' gene ', num2str(i), ' is done!']);
    end
end
end
