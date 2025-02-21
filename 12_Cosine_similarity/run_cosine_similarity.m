function run_cosine_similarity()

addpath('../generic/');


num_datasets = 8;


% get the data
if ~exist('Dall','var')
    clear Dall d
    for idataset = 1:num_datasets
        disp([num2str(idataset),'/',num2str(num_datasets)]);
        d = get_data(50,'only_tin',idataset,'dt_rel_RT',0,'positive_is_leftward',1);
        Dall(idataset) = d;
    end
end

beta_conf   = cell(num_datasets,1);
beta_choice = cell(num_datasets,1);

for i=1:num_datasets

    D = Dall(i);
    %%
    idx_neurons = 1:D.nneurons;
    S = count_spikes_in_window(D.spike_times_ms, idx_neurons, D.RT*1000-150, D.RT*1000-50)';

    % normalize??
    do_normalize = 1;
    if do_normalize
        S = nanzscore(S,[],2);
    end

    %%
    ind = D.choice==1;
    [YHAT_conf{i}, AUC_CONF(i), beta_conf{i}, AUC_CONF_stde(i)] = calc_regression_to_accuracy_crossvalid(D.correct(ind), S(:,ind));

    ind = ones(size(D.choice))==1;
    [Yhat_ch, AUC_CHOICE(i), beta_choice{i}] = calc_regression_to_accuracy_crossvalid(D.choice(ind), S(:,ind));

    YHAT_ch{i} = Yhat_ch(D.choice(ind)==1); % only the subset of Tin choices

    %%
    cosSim(i) = getCosineSimilarity(beta_conf{i}(1:end-1),beta_choice{i}(1:end-1));

    %% calc conf info along the choice axis
    ind = D.choice==1;
    [AUC_CONF_choice_axis(i), AUC_CONF_choice_axis_stde(i), AUC_boot] = ...
        auc_conf_and_bootstrap(Yhat_ch(ind), D.correct(ind));
    

end


%% calc cosyne similarity
save('./output_data/cosineSim','cosSim','beta_choice','beta_conf','do_normalize','AUC_CONF','AUC_CHOICE', ...
    'AUC_CONF_choice_axis','AUC_CONF_choice_axis_stde','AUC_CONF_stde','YHAT_conf','YHAT_ch');



end
