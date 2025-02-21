
load ../03_Kmeans_calculations/output_data/data_per_cluster.mat
%%
colores = movshon_colors(6);
p = publish_plot(1,3);
npoints = 500;
for ind=2:4
    [datamp,timeslocked] = eventlockedmatc(S(ind).H, S(ind).t, S(ind).RT ,[0.15,-0.05]);
    
    dt = S(ind).t(2)-S(ind).t(1);
    act = mean(datamp)' * 1/dt;
    
    
    p.next();
    filt = S(ind).choice==1 & S(ind).correct==1;
    plotflag = 0;
    [X,Y] = dotsanalysis.choice_against_duration_moving(act(filt),S(ind).coh(filt),S(ind).RT(filt),npoints,plotflag,ones(sum(filt),1)==1);

    for i=1:length(X)
        plot(X{i},Y{i},'Marker','none','color',colores(i,:));
        hold all
    end

    xlabel('Reaction time [s]');
    ylabel('Firing rate (sp/s)');

    title(S(ind).str);
end


p.format('FontSize',12);
same_ylim(p.h_ax);

set(gcf,'Position',[297  335  644  208]);
p.unlabel_center_plots;


p.saveas('./figures/fig_clusters_vs_coh_RT');

p.append_to_pdf('./figures/fig_clusters_vs_coh_RT');



