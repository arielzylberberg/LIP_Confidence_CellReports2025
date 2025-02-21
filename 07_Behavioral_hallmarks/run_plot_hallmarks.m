addpath('../generic/');

%%

load ../05_Regression_to_accuracy/output_data/for_hallmarks_fig.mat

%% compare with CoM data


datadir = '../data/behav_Ronald2016/';
D = load(fullfile(datadir,'pickle'));

I = D.subject==2 & D.good_trials==1;
prop_high = nanmean(D.confi(I));

pb = do_plot_conf(D.confi(I), D.coh(I), D.rt(I), D.correct(I));
set(pb.h_ax,'ylim',[0,1]);
    
%%
conf = Yhat;
p = do_plot_conf(conf, coh, RT, correct);
set(p.h_ax,'ylim',[0.6,1.001])
p.append_to_pdf('./figures/fig_hallmarks_not_thresholded',1,1);

high_conf = conf>prctile(conf,100*(1 - prop_high));
pn = do_plot_conf(high_conf, coh, RT, correct);
set(pn.h_ax,'ylim',[0,1]);


%% both in one figure

p = publish_plot(2,3);
set(gcf,'Position',[347  421  695  394]);

for i=1:3
    p.copy_from_ax(pb.h_ax(i), i);
    p.copy_from_ax(pn.h_ax(i), i+3);
end

for i=1:3
    same_xlim(p.h_ax([i,i+3]));
end

set(p.h_ax,'tickdir','out','ticklength',[0.02,0.02]);

p.ylabel('P(high confidence)',[1,4]);

% set(p.h_ax([2,3,5,6]),'xlim',[0.3,1.7]);
p.format('FontSize',11,'MarkerSize',6);

p.current_ax(3);
hl = legend_n(unique(abs(coh))*100);
set(hl,'position',[0.8498    0.7508    0.0871    0.1910]);

p.append_to_pdf('./figures/fig_hallmarks',1,1);
p.saveas('./figures/fig_hallmarks');



