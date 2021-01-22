function csn = twogenecorfuzzy(gene1, gene2, c_gene1, c_gene2, m_gene1, m_gene2, cutoff_alpha, select_id)
n2 = length(gene1);
temp = [gene1; gene2]; 
id1 = unique(c_gene1(select_id)); id2 = unique(c_gene2(select_id));
%[t_i, t_j] = find(crosstab(c_gene1, c_gene2));
d_gene1 =zeros(1, length(id2)); d_gene2 = zeros(1, length(id1));
for i = 1:length(id1)
    d_gene2(i) = std(gene2(c_gene1 == id1(i)));
end
for i = 1:length(id2)
    d_gene1(i) = std(gene1(c_gene2 == id2(i)));
end
match_id1 = find(ismember(m_gene1(1, :), id1));
match_id2 = find(ismember(m_gene2(1, :), id2));
c_gene = [repmat(m_gene1(1,match_id1), 1, length(id2)); repelem(m_gene2(1, match_id2), length(id1))];
m_gene = [repmat(m_gene1(2,match_id1), 1, length(id2)); repelem(m_gene2(2, match_id2), length(id1))];
d_gene = [repelem(d_gene1, length(id1)); repmat(d_gene2, 1, length(id2))];

upper = m_gene + d_gene;
lower = m_gene - d_gene;
csn = nan(1, n2);
for k = 1:length(id1)*length(id2)
    temp_up = repmat(upper(:, k), 1, n2);
    temp_low = repmat(lower(:, k), 1, n2);
    B = (temp - temp_up <= 0)&(temp - temp_low >= 0);
    a = sum(B,2);
    d = ((B + 0)*(B + 0)'*n2-a*a')./sqrt((a*a').*((n2-a)*(n2-a)')/(n2-1) + eps);
    csn(c_gene1 == c_gene(1, k) & c_gene2 == c_gene(2, k)) = d(1,2);
end
csn(gene1.*gene2 == 0) = -100;
if cutoff_alpha > 0
    csn = csn > -icdf('norm', cutoff_alpha, 0, 1);
end
end