function p = plot_contra_vs_ipsi()

% addpath('../functions/shadlen_script/')

%%

datadir = '../05_Regression_to_accuracy/output_data/';
a0 = load(fullfile(datadir,'auc_contra0.mat'));
a1 = load(fullfile(datadir,'auc_contra1.mat'));



V = {'TinRateRT'};
% I = find(ismember(a0.str','Tin at RT'));

p = publish_plot(1,length(V));
% set(gcf,'Position',[311  392  841  335]);


for k=1:length(V)
    p.next();
    I = find(ismember(a0.str',V{k}));


    x = a0.auc_conf(:,I);
    xse = a0.auc_conf_stde(:,I);
    
    y = a1.auc_conf(:,I);
    yse = a1.auc_conf_stde(:,I);

%     [h,pval] = ttest(x,y,'tail','left');

    [h,pval] = ttest(log(x./(1-x)),log(y./(1-y)),'tail','left');


    %%

    xlab = 'AUC_{conf}^{ipsi}';
    ylab = 'AUC_{conf}^{contra}';

    limi = [0.45,0.8];


    % error bars
    for i=1:length(x)
        plot([x(i)-xse(i), x(i)+xse(i)], [y(i),y(i)],'k');
        hold all
        plot([x(i),x(i)],[y(i)-yse(i), y(i)+yse(i)],'k');
    end
    plot(x,y,'marker','o','MarkerSize',9,'LineStyle','none','markerEdgeColor',0.0*[1 1 1],...
        'markerFaceColor',.7*[1 1 1],'LineWidth',0.1);
    

    xlim(limi)
    ylim(limi)

    axis square

    hold on

    plot(xlim,ylim,'color','k','LineStyle','--');
    xlabel(xlab)
    ylabel(ylab)
    set(gca,'tickdir','out')

    same_xytick(p.h_ax);

    ht(k) = p.text_draw(k,['p = ',num2str(redondear(pval,3))]);

end

% nontouching_spines(newScat);
p.format('FontSize',18,'LineWidthPlot',1,'LineWidthAxes',1,'MarkerSize',11);
% all_text = findall(gcf,'Type','text');
% set(all_text,'FontSize',18);


p.current_ax(1);
% hti(1) = title('With individual T_{in}^{contra} neurons');



set(ht,'FontSize',15);
displace_ax(ht,0.05,1);

% set(hti,'FontWeight','Normal','FontSize',15);

set(p.h_fig,'color','w');
set(p.h_ax,'box','off');

set(p.h_ax,'xtick',0:0.1:1);
set(p.h_ax,'ytick',0:0.1:1);

h = p.letter_the_plots();
displace_ax(h,-0.05,1);
set(h,'FontSize',20);



end
