
load ../03_Kmeans_calculations/output_data/data_per_cluster.mat
%%

K = S(1).choice==1;
n = 6;
prc=  linspace(0,100,n);

dataset = S(1).dataset(K);

aux = jitter(S(2).spk_RT(K));
Y1 = zscore_bygroup(aux, dataset);

aux = jitter(S(4).spk_RT(K));
Y3 = zscore_bygroup(aux, dataset);

[~,idx1] = index_prctile(Y1,prc);
[~,idx3] = index_prctile(Y3,prc);

[~,v1] = curva_media(Y1,idx1,[],0);
[~,v3] = curva_media(Y3,idx3,[],0);

correct = S(1).correct(K);

C = nan(n-1);
for i=1:n-1
    for j=1:n-1
        I = idx1==i & idx3==j;
        C(i,j) = nanmean(correct(I));
    end
end


p = publish_plot(1,1);
imagesc(C);
colorbar
axis xy
axis square

CT = cbrewer('div', 'RdBu', 64);
CT = CT(end:-1:1,:);
colormap(CT);

p.format();

%% smoothed
% Example matrix
data = C;

% Define original grid
[x, y] = meshgrid(1:size(data, 2), 1:size(data, 1));

% Define new high-resolution grid
[xq, yq] = meshgrid(linspace(1, size(data, 2), 100), linspace(1, size(data, 1), 100));

% Interpolate data
data_smooth = interp2(x, y, data, xq, yq, 'cubic');

% Plot the result
p = publish_plot(1,1);
set(gcf,'Position',[922  632  480  284]);
imagesc(v3,v1,data_smooth);
axis equal tight;
colormap(CT);
colorbar;
axis xy
axis square
xlabel('Firing rate, Cluster 3 (norm.)');
ylabel('Firing rate, cluster 1 (norm.)');
p.format('FontSize',12);

saveas(p.h_fig, './figures/fig_2D_clust13');

