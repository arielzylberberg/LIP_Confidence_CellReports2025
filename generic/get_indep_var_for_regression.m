function S = get_indep_var_for_regression(D, label)
% function S = get_indep_var_for_regression(D, label)

% all labels:
% labels = {
%     'sumTin',
%     'TinRateRT',
%     'TinVarRT',
%     'AllRateRT',
%     'TinToutRateRT',
%     'TinToutMinRateRT',
%     'TinMeanRT',
%     'TinMinusToutMeanRT',
%     'Tin2RateRT',
%     'Tin2MeanRT',
%     'Tout2RateRT',
%     'TinCluster1RT',
%     'TinCluster2RT',
%     'TinCluster3RT',
%     'TinCluster13RT',
%     'allButTinRateRT',
%     'TinCluster13AverageRT'
% };


do_normalize = 1;

nneurons = D.nneurons;

switch label
    case 'sumTin' % sum of Tin activity
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, 0.2, D.RT*1000-50)';

    case 'meanTin' % sum of Tin activity and average across neurons
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, 0.2, D.RT*1000-50)';
        S = nanmean(S);

    case 'TinRateRT' % firing rate at RT
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinVarRT' % variance at RT
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = nanvar(S);

    case 'AllRateRT' % firing rate at RT, all neurons
        idx_neurons = 1:nneurons;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinToutRateRT' % firing rate, Tin and Tout
        idx_neurons = [D.unitIdxLIP_Tin, D.unitIdxLIP_Tout];
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinToutMinRateRT' % firing rate, Tin, Tout, Mins
        idx_neurons = [D.unitIdxLIP_Tin, D.unitIdxLIP_Tout, D.minCells.DinRFc, D.minCells.DinRFi];
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinMeanRT' % average of Tin neurons at RT
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = nanmean(S);

    case 'TinMinusToutMeanRT' % average of the Tin minus average of the Tout
        S1 = count_spikes_in_window(D.spike_times_ms, D.unitIdxLIP_Tin, D.RT*1000-150, D.RT*1000-50)';
        S2 = count_spikes_in_window(D.spike_times_ms, D.unitIdxLIP_Tout, D.RT*1000-150, D.RT*1000-50)';
        S = nanmean(S1) - nanmean(S2);

    case 'Tin2RateRT' % firing rate at RT
        idx_neurons = D.idx_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'Tin2MeanRT' % average of Tin neurons at RT
        idx_neurons = D.idx_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = nanmean(S);

    case 'Tout2RateRT' % Tout only
        idx_neurons = D.idx_Tout;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'allButTinRateRT' % all but Tin
        ind = ~ismember(1:nneurons, D.unitIdxLIP_Tin);
        idx_neurons = find(ind);
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinCluster1RT' % firing rate at RT, subset kmeans
        clust = load('../03_Kmeans_calculations/output_data/out_3clusters.mat');
        I = ismember(clust.clusters(clust.data_id==D.idataset),[1]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinCluster2RT' % firing rate at RT, subset kmeans
        clust = load('../03_Kmeans_calculations/output_data/out_3clusters.mat');
        I = ismember(clust.clusters(clust.data_id==D.idataset),[2]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinCluster3RT' % firing rate at RT, subset kmeans
        clust = load('../03_Kmeans_calculations/output_data/out_3clusters.mat');
        I = ismember(clust.clusters(clust.data_id==D.idataset),[3]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinCluster13RT' % firing rate at RT, subset kmeans
        clust = load('../03_Kmeans_calculations/output_data/out_3clusters.mat');
        I = ismember(clust.clusters(clust.data_id==D.idataset),[1,3]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    case 'TinCluster13AverageRT' % kmeans 13, but averaging across neurons within each cluster
        clust = load('../03_Kmeans_calculations/output_data/out_3clusters.mat');
        I = ismember(clust.clusters(clust.data_id==D.idataset),[1]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        Saux = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = nanmean(Saux,1);
        I = ismember(clust.clusters(clust.data_id==D.idataset),[3]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        Saux = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = [S; nanmean(Saux,1)];

    case 'TinCluster2AverageRT' % kmeans 2, but averaging across neurons in cluster
        clust = load('../03_Kmeans_calculations/output_data/out_3clusters.mat');
        I = ismember(clust.clusters(clust.data_id==D.idataset),[2]);
        idx_neurons = D.unitIdxLIP_Tin(I);
        Saux = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = nanmean(Saux,1);

    case 'TinRateRT_odd' % firing rate at RT, odd trials
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = S(:,1:2:end);

    case 'TinRateRT_even' % firing rate at RT, even trials
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        S = S(:,2:2:end);

    case 'TinBaseline'
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, 0, 200)';

    case 'TinResidualsCohRT' % for reviewers
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

        Sr = nan(size(S));

        % calc residuals
        prc= [0:10:100];
        ucoh = nanunique(D.coh);
        [~, idxRT] = index_prctile(D.RT,prc);
        uchoice = nanunique(D.choice);
        for k=1:length(uchoice)
            for i=1:length(prc)-1
                for j=1:length(ucoh)
                    I = idxRT==i & D.coh==ucoh(j) & D.choice==uchoice(k);
                    Sr(:,I) = S(:,I) - nanmean(S(:,I),2);
                end
            end
        end
        S = Sr;

    case 'EyeBeforeAfter'
        S = [D.eyeBeforeSacc, D.eyeAfterSacc]';
        do_normalize = 0;

    case 'TinRateRT_shorter_window' % firing rate at RT
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-100, D.RT*1000-50)';

    case 'TinRateRT_plus_acoh_RT'
        idx_neurons = D.unitIdxLIP_Tin;
        S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';
        % add extra regressors: abs(coherence) and RT
        S = [S; D.RT(:)'; abs(D.coh(:))'];

end


if do_normalize
    S = nanzscore(S,[],2);
end

end
