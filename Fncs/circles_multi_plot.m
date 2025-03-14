function circles_multi_plot(Ac,As,titlename,xname,yname,xlim,ylim,legendArrary,FontSize)

% figure setting
axis equal;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xticks([ 0,Ac-As,Ac,Ac+As]);
xlabel(xname);
ylabel(yname);
title(titlename);
box on;
set(gca, 'FontName', 'Arial', 'FontSize', FontSize);
% set(gcf,'Position', [0, 0, 480, 480]);
set(gca, 'LineWidth', 1.25);
set(gca, 'YLim',ylim);
% 图示
    lgd = legend(legendArrary,'Location','best','Orientation', 'horizontal');
    set(lgd,'FontName','Arial','FontSize',FontSize);

        set(lgd, 'Color', 'none'); % 设置图例框的颜色为'none'
        set(lgd, 'Box', 'off');


end