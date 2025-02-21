addpath('../generic');

%%

% conf = load('../06_regression_to_accuracy/auc_contra1.mat');
conf = load('../05_Regression_to_accuracy/output_data/auc_contra1.mat');


clear L



L{1} = {'TinCluster1RT','TinCluster2RT','TinCluster3RT'};
L{2} = {'TinCluster13AverageRT','TinRateRT','TinMeanRT'};

% per session and average
N = length(L);
p = publish_plot(1,N);
set(gcf,'Position',[290  395  364  302]);
p.shrink(1:N,0.8,0.6);
p.displace_ax(1:N,0.15,2);
p.displace_ax(1:N,0.15,1);


for i=1:length(L)

    p.next();
    [ind, xc, mc, stdec, labels] = plot_auc_values(conf.auc_conf_all_ses, conf.auc_conf_all_ses_stde, L{i}, conf.str);


end
p.current_ax(1);
ylabel('AUC confidence');

same_xscale(p.h_ax);
set(p.h_ax(2:N),'ycolor','none');

% p.letter_the_plots();
set(p.h_ax,'ylim',[0.55,0.75]);
p.format('FontSize',14,'MarkerSize',9,'LineWidthPlot',1,'LineWidthAxes',1);
set(gca,'ytick',0.55:0.05:0.75);
drawnow
xtickangle(45);


set(p.h_ax(1),'xticklabel',{'Cluster 1','Cluster 2','Cluster 3'});
set(p.h_ax(2),'xticklabel',{'\langle Clust. 1&3 \rangle','TinCon','\langle TinCon \rangle'})


for i=1:2
    bracket_on_top_fig(p.h_ax(i), [1,1.95], 1,'add');
    bracket_on_top_fig(p.h_ax(i), [2.05,3], 1,'add');
    bracket_on_top_fig(p.h_ax(i), [1,3], 2.5,'add');
end

p.append_to_pdf('./figures/fig_auc_per_cluster',1,1);






