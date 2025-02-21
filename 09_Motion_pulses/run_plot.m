addpath('../functions/');

%%
load ../03_Kmeans_calculations/output_data/out_3clusters.mat

load('./output_data/beta_pulses');

%%
[~,ind] = sort(all_beta(:,1));

%%
clear h hfit

colores = movshon_colors(3);
p = publish_plot(3,1);
set(gcf,'Position',[743  417  334  774]);

for i=1:3
    I = clusters==i;
    p.next();

    y = nanmean(B(:,I),2);

    % normalize the weights 
    tind = ts>=0 & ts<=0.2;
    y = y - nanmean(y(tind));


    e = stderror(B(:,I),2);
    [errorPatch,dataLine] = niceBars2(ts,y,e,colores(i,:),0.3);
    delete(dataLine);
    hold all
    [c2,gof,vals(i),fallo,fit_params(i,:)] = fit_PRR_function(ts,y);

    hfit(i) = plot(vals(i).t, vals(i).fit,'color',colores(i,:));


    if i==2
        ylabel('Effect of pulse on activity (\beta_{norm})');
    end
    if i==3
        xlabel('Time from pulse onset [s]');
    end
    h(i) = plot(xlim,[0,0],'k');

end

p.format();

set(p.h_ax,'xlim',[0,0.8]);

% p.unlabel_center_plots();
p.format('LineWidthAxes',1,'FontSize',20);
set(hfit,'LineWidth',2);
set(h,'LineWidth',1);

p.append_to_pdf('./figures/fig_pulses_per_clust_baseline');

save('./output_data/out_fit','fit_params');

