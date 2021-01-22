function [] = temp_avg_heat(cdata1, cdata2, xvalues, yvalues, color_str, title_str, sup_title)
figure;
subplot(1, 2, 1); [~, n] = size(cdata1);
for i = 1:n
    cdata1(i, i) = NaN;
end
h = heatmap(xvalues,yvalues, cdata1);
colormap(brewermap([],color_str));
colorbar
h.Title = [title_str ', Group 1'];
h.XLabel = 'Genes';
h.YLabel = 'Genes';
subplot(1, 2, 2);
for i = 1:n
    cdata2(i, i) = NaN;
end
h = heatmap(xvalues,yvalues, cdata2);
colormap(brewermap([],color_str));
colorbar
h.Title = [title_str ', Group 2'];
h.XLabel = 'Genes';
h.YLabel = 'Genes';
sgtitle(sup_title);
end