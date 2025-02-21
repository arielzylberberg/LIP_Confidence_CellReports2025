function [Yhat, AUC, beta, AUC_stde, pval, Yhat_not_crossvalid, AUC_not_crossvalid, beta_se] = ...
    calc_regression_to_accuracy_crossvalid(correct, neural_activity, do_boot)

if nargin<3 || isempty(do_boot)
    do_boot = 1;
end


nfolds = 10;
N = length(correct);
c = cvpartition(N,'KFold',nfolds);
% acc_pred = 0;
Yhat = nan(N,1);
for i=1:nfolds

    % fit
    I = c.training(i)==1;
    indepvar = {'x',neural_activity(:,I)','ones',ones(sum(I),1)};
    depvar = correct(I);
    [beta,idx,stats,x,LRT] = f_regression(depvar,[],indepvar);

    % test
    I = c.test(i)==1;
    x_test = [neural_activity(:,I)',ones(sum(I),1)];
    yhat = glmval(beta,x_test,'logit','constant','off');

%     acc_pred = acc_pred + mean((yhat>0.5)==correct(I));
    Yhat(I) = yhat;

    Betas(i,:) = beta;
end



%%
SCORES = Yhat;
LABELS = correct;
POSCLASS = 1;
[X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(LABELS,SCORES,POSCLASS);

%% bootrstrap the auc
if do_boot
    nboot = 5000;
    AUC_boot = nan(nboot,1);
    n = length(correct);
    I = ceil(rand(n,nboot)*n);
    for i=1:nboot
        SCORES = Yhat(I(:,i));
        LABELS = correct(I(:,i));
        [~,~,~,AUC_boot(i)] = perfcurve(LABELS,SCORES,POSCLASS);
    end
    AUC_stde = std(AUC_boot);
else
    AUC_stde = [];
end

%%
% beta = nanmean(Betas); % average over the folds

%% do a fit with all trials for the beta and pval
I = ones(size(correct))==1;
indepvar = {'x',neural_activity(:,I)','ones',ones(sum(I),1)};
depvar = correct(I);
[beta,~,stats,x] = f_regression(depvar,[],indepvar);
pval = stats.p;
Yhat_not_crossvalid = glmval(beta,x,'logit','constant','off');

beta_se = stats.se;

SCORES = Yhat_not_crossvalid;
LABELS = correct;
POSCLASS = 1;
[~,~,~,AUC_not_crossvalid] = perfcurve(LABELS,SCORES,POSCLASS);


end
