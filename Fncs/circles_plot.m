function circles_plot(Ac,As,titlename,xname,yname,Xlim,Ylim,legendArrary,FontSize)
viscircles([Ac 0],As,'Color','black',LineWidth=1.5);
% figure setting
axis equal;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% 绘制x，y轴 带箭头的 
pos = ax.Position;

% 计算绘图区域的宽度和高度
width = pos(3);
height = pos(4);

% 计算箭头的起点和终点坐标
x_start = pos(1);
x_end = pos(1) + width;
y_start = pos(2);
y_end = pos(2) + height;
% 计算y轴中点
y=(y_start+y_end)/2;
% 绘制x轴箭头
annotation('arrow', [x_start, x_end], [y, y], 'LineWidth', 1.5);
% 绘制y轴箭头
annotation('arrow', [x_start, x_start], [y_start, y_end], 'LineWidth', 1.5);


xticks([ 0,Ac-As,Ac,Ac+As]);

title(titlename);

set(gca, 'FontName', 'Arial', 'FontSize', FontSize);
set(gca, 'FontName', '宋体');
% set(gcf,'Position', [0, 0, 480, 480]);
set(gca, 'LineWidth', 1.25);
% set(gca,'XLim',Xlim);
set(gca, 'YLim',Ylim);
% 隐藏默认的坐标轴和刻度
set(gca, 'YTick', []);

xticks(4);
yticks(0);


xlabel(xname);
ylabel(yname,'Rotation',0,'Position',[x_end+0.3, y_end+0.7]);
% 图示
    lgd = legend(legendArrary,'Location','best');
    set(lgd,'FontName','Arial','FontSize',FontSize);

        set(lgd, 'Color', 'none'); % 设置图例框的颜色为'none'
        set(lgd, 'Box', 'off');


end