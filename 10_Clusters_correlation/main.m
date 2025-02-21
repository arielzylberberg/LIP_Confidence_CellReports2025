addpath('../generic');
addpath('./functions'); 


%%
do_save_flag = 1;

win = 0.025;
S = get_and_prep_data(win);


color_scale = 3; % 4
filt_flag = 1; % 1

%%
xind = [3,4,4];
yind = [2,2,3];
XLABEL = {'Time [s], cluster 2',...
    'Time [s], cluster 3',...
    'Time [s], cluster 3',...
    };
YLABEL = {'Time [s], cluster 1',...
    'Time [s], cluster 1',...
    'Time [s], cluster 2',...
    };


RT = S(1).RT;
coh = S(1).coh;
choice = S(1).choice;
t = S(1).t;

p = publish_plot(3,1);
set(gcf,'Position',[982   345   276   790]);


plot_order = [1:3];
for k = 1:length(xind)

    X = S(xind(k)).zHres;
    Y = S(yind(k)).zHres;
    
    
    switch color_scale
        case 1
            colores = cbrewer('div','RdBu',100);
        case 2
            colores = cbrewer('seq','YlGnBu',100);
            colores = colores(end:-1:1,:);
        case 3
            colores = cbrewer('div','RdBu',100);
            colores = colores(end:-1:1,:);
        case 4
            colores = cbrewer('div','RdBu',100);
            colores = colores(end:-1:1,:);

    end
    

    switch filt_flag
        case 1
%             I = choice==0;
            I = choice==1;
            tind = find(t>=0.2 & t<=0.8);
        case 2
            I = RT>0.55;
            tind = find(t>=0.0 & t<=0.5);
    end


    if 0
        nt = length(tind);
        rho = nan(nt);
        pval = nan(nt);
        for i=1:nt
            for j=1:nt
                x = X(:,tind(i));
                y = Y(:,tind(j));
                K = ~isnan(x) & ~isnan(y) & I;
                [rho(i,j),pval(i,j)] = corr(x(K),y(K));
            end
        end
    end

    [rho, pval, p_as_extreme(k)] = do_perm_test_corr(X(I,:),Y(I,:),tind);


    p.current_ax(plot_order(k));

    switch color_scale
        case 1

            lim = max(abs(rho(:)));
            lim = [-lim,lim];
            imagesc(t(tind),t(tind),rho',lim);
        case 2
            imagesc(t(tind),t(tind),rho');
        case 3
            lim = max(abs(rho(:)));
            lim = [-lim,lim];
            imagesc(t(tind),t(tind),rho',lim);
        case 4
            lim = [-0.45,0.45];
            imagesc(t(tind),t(tind),rho',lim);

    end
    axis square

    colormap(colores);
    colorbar 

    axis xy
    h = refline(1,0);
    set(h,'color','k','LineStyle','--');
    xlabel(XLABEL{k});
    ylabel(YLABEL{k});


end

p.format('FontSize',12);
p.append_to_pdf('./figures/fig_corr_after_grouping_across_sessions',1,do_save_flag);
saveas(p.h_fig,'./figures/fig_corr_after_grouping_across_sessions');

save('./output_data/stats_corr','p_as_extreme','XLABEL','YLABEL');

%%

