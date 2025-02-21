addpath('../generic');
addpath('functions');

%%

p1 = plot_contra_vs_ipsi();
p1.append_to_pdf('./figures/fig_contra_vs_ipsi',1,1);

p2 = illustrate_auc_method();
p2.append_to_pdf('./figures/fig_illustrate_AUC',1,1);

p = publish_plot(1,2);
set(gcf,'Position',[412  605  718  297]);
p.new_axes('Position',[0.15,0.5,0.15,0.4]);


p.copy_from_ax(p1.h_ax(1),2);
p.copy_from_ax(p2.h_ax(1),1);

p.copy_from_ax(p2.h_ax(2),3);

p.current_ax(2);
axis square

p.current_ax(3);
axis square



set(p.h_ax(1),'ycolor','none');

p.format('MarkerSize',[10,10],'LineWidthAxes',1,'FontSize',16);

set(p.h_ax(3),"FontSize",11);

set(p.h_ax(2),'xtick',0.5:0.1:0.8,'ytick',0.5:0.1:0.8);

h = p.letter_the_plots('show',1:2);
set(h,'FontSize',17);

p.append_to_pdf('./figures/fig_combined',1,1);