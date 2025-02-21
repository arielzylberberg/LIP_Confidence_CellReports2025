function legend_dots_coh(colores)

figure();
delta = 0.1;
str = {'0','3.2','6.4','12.8','25.6','51.2'};
for i=1:6
    h = rectangle('Position',[0+delta,i+delta,1-delta*2,1-delta*2],'FaceColor',colores(7-i,:),'LineStyle','none');
    hold all

    h = rectangle('Position',[1.5+delta,i+delta,1-delta*2,1-delta*2],'FaceColor',colores(i+6,:),'LineStyle','none');

    hf(i) = text(3, i+0.5, str{i},'horizontalalignment','left','verticalalignment','middle');
end

text(2, 8, 'Motion strength [%coh]','horizontalalignment','center','verticalalignment','middle','FontSize',20);

xlim([-1,10])
ylim([-1,10])
axis square
axis off
set(hf,'FontSize',20);

end