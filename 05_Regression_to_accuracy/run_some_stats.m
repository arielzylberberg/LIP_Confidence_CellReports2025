load ./auc_contra1.mat

%%

ind1 = ismember(str,'TinCluster1RT');
ind2 = ismember(str,'TinCluster2RT');
ind3 = ismember(str,'TinCluster3RT');

nboot = size(AUC_boot_all_ses);
idx = combvec(1:nboot,1:nboot);

%% stats comparison clusters
m_12 = mean(AUC_boot_all_ses(idx(:,1),ind1)>AUC_boot_all_ses(idx(:,2),ind2));
m_32 = mean(AUC_boot_all_ses(idx(:,1),ind3)>AUC_boot_all_ses(idx(:,2),ind2));

%% another comparison

ind1 = ismember(str,'TinRateRT');
ind2 = ismember(str,'TinCluster13AverageRT');

nboot = size(AUC_boot_all_ses);
idx = combvec(1:nboot,1:nboot);

m_12 = mean(AUC_boot_all_ses(idx(:,1),ind1)>AUC_boot_all_ses(idx(:,2),ind2));

%% eye
ind1 = ismember(str,'TinRateRT');
ind2 = ismember(str,'EyeBeforeAfter');

nboot = size(AUC_boot_all_ses);
idx = combvec(1:nboot,1:nboot);

m = mean(AUC_boot_all_ses(idx(:,1),ind1)>AUC_boot_all_ses(idx(:,2),ind2));


