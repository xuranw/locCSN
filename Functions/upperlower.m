function [upper, lower] = upperlower(data, boxsize)
[n1,n2] = size(data);
upper = zeros(n1,n2);
lower = zeros(n1,n2);
for i = 1:n1
    [s1,s2] = sort(data(i,:));
%    n0 = n2-sum(sign(s1));
    h = round(boxsize/2*n2);
    k = 1;
    while k <= n2
        s = 0;
        while k+s+1 <= n2 && s1(k+s+1) == s1(k)
            s = s+1;
        end
        if s >= h
            upper(i,s2(k:k+s)) = data(i,s2(k));
            lower(i,s2(k:k+s)) = data(i,s2(k));
        else
            upper(i,s2(k:k+s)) = data(i,s2(min(n2,k+s+h)));
            lower(i,s2(k:k+s)) = data(i,s2(max(1,k-h)));
        end
        k = k+s+1;
    end
end
end