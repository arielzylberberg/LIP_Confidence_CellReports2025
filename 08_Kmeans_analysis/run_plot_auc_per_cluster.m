addpath('../generic');

%%

% conf = load('../0    _regression_to_accuracy/auc_contra1.mat');
conf = load('../05_Regression_to_accuracy/output_data/auc_contra1.mat');
clear L


L =  {'TinCluster1RT','TinCluster2RT','TinCluster3RT',};


% per session and average
p = publish_plot(1,1);
set(gcf,'Position',[1060   786   224   341]);

p.next();
[ind, xc, mc, stdec, labels] = plot_auc_values(conf.auc_conf_all_ses, conf.auc_conf_all_ses_stde, L, conf.str);
ylabel('AUC confidence');

ylim([0.55,0.75]);
p.format('FontSize',14,'MarkerSize',10,'LineWidthPlot',1,'LineWidthAxes',0.5);
set(gca,'ytick',0.55:0.05:0.75);
xtickangle(35);
drawnow



p.append_to_pdf('./figures/fig_auc_per_cluster',1,1);


