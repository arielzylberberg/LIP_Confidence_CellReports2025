function p = change_colors_PSTH(p)

for j=1:length(p.h_ax)
    h  = get(p.h_ax(j),'children');


    colores = movshon_colors(6);
    %             colores = [colores;colores(end:-1:1,:)];
    colores = [colores(end:-1:1,:); colores];
    lsty = {'-',':'};


    for i=1:length(h)
        set(h(i),'color',colores(i,:));
        if i<=6
            set(h(i),'LineStyle',lsty{1});
        else
            set(h(i),'LineStyle',lsty{2});
        end
    end
end


for i=1:2:length(p.h_ax)
    p.current_ax(i);
    hold all
    plot([0,0],ylim,'color',0.7*[1,1,1]);
end


set(p.h_ax(1:2:end),'xtick',[0,0.25,0.5,0.75],'xticklabel',[0,0.25,0.5,0.75]);
set(p.h_ax(2:2:end),'xtick',[-0.5:0.25:0],'xticklabel',[-0.5:0.25:0],'ycolor','none');
set(p.h_ax,'tickdir','out','ticklength',[0.02,0.02]);
p.unlabel_center_plots();
p.format('FontSize',12);

end