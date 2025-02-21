pd = publish_plot(4,2);
pd.load_from_fig('../03_Kmeans_calculations/figures/fig_kmeans.fig');
set(pd.h_fig,'Position',[440  180  400  618]);

for j=1:length(pd.h_ax)
    h  = get(pd.h_ax(j),'children');

    color_flag = 1;
    switch color_flag
        case 1
            colores = movshon_colors(6);
%             colores = [colores;colores(end:-1:1,:)];
            colores = [colores(end:-1:1,:); colores];
            lsty = {'-',':'};
        case 2
            colores = cbrewer('div','BrBG',12);
            lsty = {'-','-'};
        case 3
            colores = cbrewer('div','RdYlBu',12);
            lsty = {'-','-'};
        case 4
            colores = NS_colors(12);
            lsty = {'-','-'};

    end

    for i=1:length(h)
        set(h(i),'color',colores(i,:));
        if i<=6
            set(h(i),'LineStyle',lsty{1});
        else
            set(h(i),'LineStyle',lsty{2});
        end
    end
end

set(pd.h_ax,'ylim',[5,30]);
set(pd.h_ax(3:4),'ylim',[5,32]);
set(pd.h_ax(7:8),'ylim',[5,50]);

for i=1:2:length(pd.h_ax)
    pd.current_ax(i);
    hold all
    plot([0,0],ylim,'color',0.7*[1,1,1]);
end


set(pd.h_ax(1:2:end),'xtick',[0,0.25,0.5,0.75],'xticklabel',[0,0.25,0.5,0.75]);
set(pd.h_ax(2:2:end),'xtick',[-0.5:0.25:0],'xticklabel',[-0.5:0.25:0],'ycolor','none');
set(pd.h_ax,'tickdir','out','ticklength',[0.02,0.02]);
pd.unlabel_center_plots();
pd.format('FontSize',12');

saveas(pd.h_fig, './figures/fig_kmeans_recolored.fig');
