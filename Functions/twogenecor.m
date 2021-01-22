function csn = twogenecor(gene1, gene2, c_gene1, c_gene2, select_id, cutoff_alpha, fuzzy, m_gene1, m_gene2)
if nargin < 7 || isempty(m_gene1) || isempty(m_gene2)
  fuzzy = false;
end
if nargin <6 || isempty(cutoff_alpha)
   cutoff_alpha = 0.01;
end
if nargin < 5 || isempty(select_id)
    select_id = 1:length(gene1);
end
if fuzzy
    csn = twogenecorfuzzy(gene1, gene2, c_gene1, c_gene2, m_gene1, m_gene2, cutoff_alpha, select_id);
else
    n2 = length(gene1);
    temp = [gene1; gene2]; 
    d_gene1 = zeros(1, n2); d_gene2 = zeros(1, n2);
    id1 = unique(c_gene1); id2 = unique(c_gene2);
    %[t_i, t_j] = find(crosstab(c_gene1, c_gene2));
    for i = id1
        d_gene2(c_gene1 == i) = std(gene2(c_gene1 == i));
    end
    for i = id2
        d_gene1(c_gene2 == i) = std(gene1(c_gene2 == i));
    end
    upper = [gene1 + d_gene1; gene2 + d_gene2];
    lower = [gene1 - d_gene1; gene2 - d_gene2];
    csn = nan(1, n2);
    for k = select_id'
        if gene1(k)*gene2(k) == 0
            csn(k) = -100;
        else
            temp_up = repmat(upper(:, k), 1, n2);
            temp_low = repmat(lower(:, k), 1, n2);
            B = (temp - temp_up <= 0)&(temp - temp_low >= 0);
            a = sum(B,2);
            d = ((B + 0)*(B + 0)'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1) + eps);
            csn(k) = d(1,2);
        end
    end
    if cutoff_alpha > 0
        csn = csn > -icdf('norm', cutoff_alpha, 0, 1);
    end
end
end