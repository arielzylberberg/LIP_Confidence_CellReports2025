function [clusters, all_beta, all_stats] = do_regre_and_kmeans(spk, n_clusters, correct, RT, choice, coh, neuron_id)

u = nanunique(neuron_id);

nneurons = length(u);

for i=1:length(u)

    %% compute residuals from regression to coherence

    I = neuron_id == u(i) & correct == 1;
    depvar = spk(I);

    % zscore
    depvar = nanzscore(depvar);


    indepvars = {'coh',coh(I),'RT',RT(I),'choice',choice(I),'ones',ones(sum(I),1)};
    [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvars);

    nb = length(indepvars)/2 - 1; % number of rel. variables in the regression (minus the offset)

    if i==1
        all_beta = nan(nneurons, nb);
        all_stats = nan(nneurons,nb);
    end

    all_beta(i,:) = beta(1:nb);
    all_stats(i,:) = stats.p(1:nb);

end


%% do the clustering

I = ones(nneurons, 1)==1;
n = size(all_beta,1);
clusters = nan(n,1);

rng(26,'twister');

clusters(I) = kmeans(all_beta(I,:),n_clusters);

%% resort the clusters by some criteria...
% criteria: average first regre coeff

figure();
[~,xx] = curva_media(all_beta(:,1), clusters,~isnan(clusters),0);

% [~,ind] = sort(xx);
[~,ind] = sort(xx,'descend');

% resorting
clusters_resort = nan(size(clusters));
for i=1:n_clusters
    I = clusters==ind(i);
    clusters_resort(I) = i;
end
clusters = clusters_resort;







