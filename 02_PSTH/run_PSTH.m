addpath('../generic/');


%%
load('../data_prepro/LIP_data_concat/data_Tin_with_H.mat','t','H','RT','coh_extra','correct');

%%
[p, av] = fn_plot_PSTH(t,H,RT,coh_extra,correct);

set(p.h_ax,'ylim',[5,30]);

p.append_to_pdf('fig_PSTH',1,1);
