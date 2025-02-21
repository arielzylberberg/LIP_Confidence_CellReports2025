function regression_to_conf_group_sessions(labels, contralateral_choice_flag)


if nargin==1 || isempty(contralateral_choice_flag)
    regression_to_conf_group_sessions(labels, 1);
    regression_to_conf_group_sessions(labels, 0);
    return;
end


nlabels = length(labels);



num_datasets = 8;
auc_conf = nan(num_datasets, nlabels);
auc_conf_stde = auc_conf;
auc_conf_all_ses = nan(nlabels,1);
auc_conf_all_ses_stde = nan(nlabels,1);

AUC_boot_all_ses = [];

% get the data
if ~exist('Dall','var')
    clear Dall d
    for idataset = 1:num_datasets
        disp([num2str(idataset),'/',num2str(num_datasets)]);
        d = get_data(50,[],idataset,'dt_rel_RT',0,'positive_is_leftward',1);
        Dall(idataset) = d;
    end
end


pred = struct();
filename = ['auc_contra',num2str(contralateral_choice_flag)];
%

for k = 1:nlabels


    v_yhat = [];
    correct = [];
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

            switch contralateral_choice_flag
                case 1
                    tr_ind = D.choice==1;

                case 0
                    tr_ind = D.choice==0;

            end
            
            %% prep
            dat.x = S(:,tr_ind);

            dat.y = D.correct(tr_ind)==1;
            dat.RT = D.RT(tr_ind);
            dat.coh = D.coh(tr_ind);
            dat.correct = D.correct(tr_ind);

            %% do regression
            [Yhat, AUC, beta, AUC_stde] = calc_regression_to_accuracy_crossvalid(dat.correct, dat.x);
            auc_conf(idataset,k) = AUC;
            auc_conf_stde(idataset,k) = AUC_stde;
            Betas{idataset,k} = beta;

            %% save session data
            v_yhat  = [v_yhat; Yhat];
            correct = [correct; dat.correct];
            coh     = [coh; dat.coh];
            RT      = [RT; dat.RT];
            NFactors = [NFactors; size(dat.x,1)];

        end

    end

    Yhat = v_yhat;

    %% auc after grouping across sessions
    SCORES = Yhat;
    LABELS = correct;


    num_factors(k) = nanmean(NFactors);

    % bootstrap to get se?
    [auc_conf_all_ses(k), auc_conf_all_ses_stde(k), AUC_boot_all_ses(:,k)] = auc_conf_and_bootstrap(SCORES, LABELS);
    

    %%

    p = do_plot_conf(Yhat, coh, RT, correct);

    p.text_draw_fig(labels{k});
    p.append_to_pdf(['./figures/fig_',filename],k==1,1);

    % save the one for the paper
    if isequal(labels{k},'TinRateRT') && contralateral_choice_flag==1
        p.saveas(['./figures/fig_conf_grouped_c',num2str(contralateral_choice_flag)]);

        % save the data to re-do the figure 
        save('./output_data/for_hallmarks_fig','Yhat','coh','RT','correct');

        % also save to illustrate AUC method
        save('./output_data/for_illustration_AUC','SCORES','LABELS');
    
    end

    pred(k).Yhat = Yhat;
    pred(k).coh = coh;
    pred(k).RT = RT;
    pred(k).correct = correct;

end

%%
str = labels;
save(fullfile('output_data',filename),'auc_conf','auc_conf_all_ses','auc_conf_all_ses_stde','str',...
    'num_factors','auc_conf_stde','AUC_boot_all_ses','pred','Betas');


%%
% per session and average
p = publish_plot(1,1);
set(gcf,'Position',[514  409  412  370]);
[~,ind] = sort(nanmean(auc_conf),'ascend');
terrorbar(1:length(ind), nanmean(auc_conf(:,ind)), stderror(auc_conf(:,ind)),'marker','o','markerfacecolor','w');
grid on
set(gca,'xticklabel',str(ind),'xtick',1:length(ind));
ylabel('AUC confidence');
p.format('FontSize',11);
p.append_to_pdf(['./figures/fig_',filename],0,1);


% doing the auc on all sessions
p = publish_plot(1,1);
set(gcf,'Position',[514  409  412  370]);
[~,ind] = sort(auc_conf_all_ses,'ascend');
plot(1:length(ind), auc_conf_all_ses(ind),'marker','o','markerfacecolor','w');
grid on
set(gca,'xticklabel',str(ind),'xtick',1:length(ind));
ylabel('AUC confidence');
p.format('FontSize',11);
p.append_to_pdf(['./figures/fig_',filename],0,1);


end

