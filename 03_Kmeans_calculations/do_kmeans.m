function do_kmeans()


addpath('../generic');

%%

rng(2424890,'twister');

n_clusters = 3; 


%%

load('../data_prepro/LIP_data_concat/data_Tin_with_H.mat','H','RT','choice','coh','coh_extra','correct','dataset','neuron_id','spk_RT','t');

%%

[clusters, all_beta, all_stats] = do_regre_and_kmeans(spk_RT, n_clusters, correct, RT, choice, coh, neuron_id);
[~,data_id] = curva_media(dataset, neuron_id,[],0);


%% plot on pca axes
% pca on the beta coeff
[COEFF, SCORE, LATENT, TSQUARED] = pca(all_beta);
u = nanunique(clusters);
p = publish_plot(1,1);
for i=1:length(u)
    I = clusters==i;

    plot(SCORE(I,1),SCORE(I,2),'.');
    hold all
    xlabel('PCA 1');
    ylabel('PCA 2');
    hl = legend_n(1:length(u));
    set(hl,'location','best');
    p.format('MarkerSize',13);

end

p.append_to_pdf('./figures/fig_pca_single',1,1);


save(['./output_data/out_',num2str(n_clusters),'clusters'],'clusters','data_id','all_beta');


%% clusters vs. betas

u = nanunique(clusters);
p = publish_plot(1,3);
lab = {'coh','RT','choice'};
comps = [1,2;2,3;3,1];
for j=1:size(comps,1)
    p.next();
    for i=1:length(u)
        I = clusters==i;
        plot(all_beta(I,comps(j,1)), all_beta(I,comps(j,2)),'.');
        hold all
        xlabel(['beta ',lab{comps(j,1)}]);
        ylabel(['beta ',lab{comps(j,2)}]);
    end
end
hl = legend_n(1:length(u));
set(hl,'location','best');
p.format('MarkerSize',13);

p.append_to_pdf('./figures/fig_cluster_by_var',1,1);


%% 3-d plot 
I = ~isnan(clusters);
colores = movshon_colors(3);
color = colores(clusters(I),:);
p = scatter3_with_projections(all_beta(I,:), color);
xlabel('\beta_{coh}');
ylabel('\beta_{RT}');
zlabel('\beta_{choice}');
p.format();
set(gcf,'Position',[584  591  408  410]);
p.saveas('./figures/fig_3d_cluster');
p.append_to_pdf('./figures/fig_3d_cluster',1,1);

%% plot the clusters

% correct trials
p = fn_plot_PSTH_per_cluster(t,H,RT,coh_extra, correct, clusters, neuron_id);
p.append_to_pdf('./figures/fig_kmeans_alltogether_single',1,1);
p.saveas('./figures/fig_kmeans');

% error trials
p = fn_plot_PSTH_per_cluster(t,H,RT,coh_extra, correct==0, clusters, neuron_id);
p.append_to_pdf('./figures/fig_kmeans_alltogether_single',0,1);

% contra vs ipsi
p = fn_plot_PSTH_per_cluster(t,H,RT,choice, ismember(correct,[0,1]), clusters, neuron_id);
p.current_ax(1);
hl = legend('contra','ipsi');
set(hl,'location','best');
p.append_to_pdf('./figures/fig_kmeans_alltogether_single',0,1);

% split by RT
RT_median = index_prctile_by_group(RT,[0,33,67,100],ones(size(coh_extra)));
p = fn_plot_PSTH_per_cluster(t,H,RT,RT_median, choice==1, clusters, neuron_id);
p.current_ax(1);
hl = legend('fast','intermediate','slow');
set(hl,'location','best');
p.append_to_pdf('./figures/fig_kmeans_alltogether_single',0,1);
p.saveas('./figures/fig_kmeans_by_RT');

% correct and incorrect, contra only (errors grouped; coh=0 plot with
% correct ones)
I = choice==1;
c = coh(I); 
c(c<0) = -1;%just an identifier
p = fn_plot_PSTH_per_cluster(t,H(I,:),RT(I),c, ones(size(correct(I))), clusters, neuron_id(I));
p.saveas('./figures/fig_kmeans_errors_contra');
p.append_to_pdf('./figures/fig_kmeans_alltogether_single',0,1);






