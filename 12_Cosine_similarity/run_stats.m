load('./output_data/cosineSim');

%%
mean(abs(cosSim))
stderror(abs(cosSim))

%%

x = AUC_CONF;
xse = AUC_CONF_stde;

y = max(AUC_CONF_choice_axis, 1-AUC_CONF_choice_axis);
yse = AUC_CONF_choice_axis_stde;

[h,pval] = ttest(log(x./(1-x)),log(y./(1-y)));

mean(y)
stderror(y)
pval