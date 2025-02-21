function [AUC, AUC_stde, AUC_boot] = auc_conf_and_bootstrap(ConfEstimate, correct)

if nargin<3
    session_id = ones(size(correct));
end

SCORES = ConfEstimate;
LABELS = correct;
POSCLASS = 1;
[X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(LABELS,SCORES,POSCLASS);

%% bootrstrap the auc. Bootstrap within session. 

nboot = 5000;
AUC_boot = nan(nboot,1);
n = length(correct);
I = ceil(rand(n,nboot)*n);
for i=1:nboot
    
    SCORES = ConfEstimate(I(:,i));
    LABELS = correct(I(:,i));

    [~,~,~,AUC_boot(i)] = perfcurve(LABELS,SCORES,POSCLASS);

end
AUC_stde = std(AUC_boot);

end