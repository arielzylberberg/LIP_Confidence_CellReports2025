load('./output_data/cosineSim');

%% cosine sim

p = publish_plot(1,2);
% set(gcf,'Position',[737  712  494  341]);
p.next();

abs_flag = 0;

if abs_flag==0
    h = histogram(cosSim,[-1:0.05:1]);
    set(h,'FaceColor',0.7*[1,1,1]);
    hold on
    plot([0,0],[0,2.2],'k--');
    xlabel('Cosine similarity');
    xlim([-1,1]);
else
    h = histogram(abs(cosSim),[0:0.025:1]);
    set(h,'FaceColor',0.7*[1,1,1]);
    xlabel('|Cosine similarity|');
    xlim([0,1]);
end

ylabel('# sessions');

p.format();
ylim([0,2.2]);

%%
p.next();

x = AUC_CONF;
xse = AUC_CONF_stde;

y = max(AUC_CONF_choice_axis, 1-AUC_CONF_choice_axis);
yse = AUC_CONF_choice_axis_stde;

[h,pval] = ttest(log(x./(1-x)),log(y./(1-y)),'tail','right');


%%

xlab = 'AUC_{conf}^{contra}';
ylab = 'AUC_{conf}^{contra} along choice direction';

limi = [0.45,0.8];


% error bars
for i=1:length(x)
    plot([x(i)-xse(i), x(i)+xse(i)], [y(i),y(i)],'k');
    hold all
    plot([x(i),x(i)],[y(i)-yse(i), y(i)+yse(i)],'k');
end
plot(x,y,'marker','o','MarkerSize',9,'LineStyle','none','markerEdgeColor',0.0*[1 1 1],...
    'markerFaceColor',0.8*[1 1 1],'LineWidth',0.1);


xlim(limi)
ylim(limi)

axis square

hold on

plot(xlim,ylim,'color','k','LineStyle','--');
xlabel(xlab)
ylabel(ylab)
set(gca,'tickdir','out')

same_xytick(p.h_ax(2));

ht = p.text_draw(2,['p = ',num2str(redondear(pval,3))]);

%%
p.format('FontSize',15,'LineWidthAxes',1,'LineWidthPlot',1);
p.append_to_pdf('./figures/fig_cosine',1,1);