function main()


rng(313913,'twister');

addpath('../generic/');


labels = {
    'TinRateRT',
    'TinMeanRT',
    'TinCluster1RT',
    'TinCluster2RT',
    'TinCluster3RT',
    'TinCluster13RT',
    'TinCluster13AverageRT'
};

nlabels = length(labels);


num_datasets = 8;
auc_choice = nan(num_datasets, nlabels);
auc_choice_stde = auc_choice;
auc_choice_all_ses = nan(nlabels,1);
auc_choice_all_ses_stde = nan(nlabels,1);

% get the data
if ~exist('Dall','var')
    clear Dall d
    for idataset = 1:num_datasets
        disp([num2str(idataset),'/',num2str(num_datasets)]);
        d = get_data(50,[],idataset,'dt_rel_RT',0,'positive_is_leftward',1);
        Dall(idataset) = d;
    end
end


filename = ['auc_choice'];
%

for k = 1:nlabels


    v_yhat = [];
    correct = [];
    choice = [];
    coh = [];
    RT = [];
    NFactors = [];

    disp(labels{k});

    for idataset = 1:num_datasets

        disp([num2str(idataset),'/',num2str(num_datasets)]);

        D = Dall(idataset);

        S = get_indep_var_for_regression(D, labels{k});


        %%
        if size(S,1)>0

            tr_ind = ismember(D.choice, [0,1]);

            %% prep
            dat.x = S(:,tr_ind);

            dat.y = D.correct(tr_ind)==1;
            dat.RT = D.RT(tr_ind);
            dat.coh = D.coh(tr_ind);
            dat.correct = D.correct(tr_ind);
            dat.choice = D.choice(tr_ind);

            %% do regression

            [Yhat, AUC, beta, AUC_stde] = calc_regression_to_accuracy_crossvalid(dat.choice, dat.x);
            auc_choice(idataset,k) = AUC;
            auc_choice_stde(idataset,k) = AUC_stde;
            Betas{idataset,k} = beta;

            %% save session data
            v_yhat  = [v_yhat; Yhat];
            correct = [correct; dat.correct];
            choice  = [choice; dat.choice];
            coh     = [coh; dat.coh];
            RT      = [RT; dat.RT];
            NFactors = [NFactors; size(dat.x,1)];

        end

    end

    Yhat = v_yhat;

    %% auc after grouping across sessions
    SCORES = Yhat;
    LABELS = choice;


    num_factors(k) = nanmean(NFactors);

    % bootstrap to get se?
    [auc_choice_all_ses(k), auc_choice_all_ses_stde(k)] = auc_conf_and_bootstrap(SCORES, LABELS);


    %%
    Yhat_choice(:,k) = v_yhat; 

end

%%
str = labels;
save(['./output_data/',filename],'auc_choice','auc_choice_all_ses','auc_choice_all_ses_stde',...
    'str','num_factors','auc_choice_stde','Yhat_choice','Betas');


%%
% per session and average
p = publish_plot(1,1);
set(gcf,'Position',[514  409  412  370]);
[~,ind] = sort(nanmean(auc_choice),'ascend');
terrorbar(1:length(ind), nanmean(auc_choice(:,ind)), stderror(auc_choice(:,ind)),'marker','o','markerfacecolor','w');
grid on
set(gca,'xticklabel',str(ind),'xtick',1:length(ind));
ylabel('AUC choice');
p.format('FontSize',11);
p.append_to_pdf(['./figures/fig_',filename],0,1);


% doing the auc on all sessions
p = publish_plot(1,1);
set(gcf,'Position',[514  409  412  370]);
[~,ind] = sort(auc_choice_all_ses,'ascend');
plot(1:length(ind), auc_choice_all_ses(ind),'marker','o','markerfacecolor','w');
grid on
set(gca,'xticklabel',str(ind),'xtick',1:length(ind));
ylabel('AUC choice');
p.format('FontSize',11);
p.append_to_pdf(['./figures/fig_',filename],0,1);

end

