function p = illustrate_auc_method()


load ../05_Regression_to_accuracy/output_data/for_illustration_AUC.mat

%%

I = LABELS==1;
x = SCORES(I);

I = LABELS==0;
y = SCORES(I);

p = publish_plot(1,2);
p.shrink(2,0.8,0.8);

p.next();
h1 = histogram(x);
hold on
h2 = histogram(y);

h1.Normalization = 'probability';
h1.BinWidth = 0.02;
h2.Normalization = 'probability';
h2.BinWidth = 0.02;

set([h1,h2],'EdgeColor','none');
hl = legend('Correct','Incorrect');
set(hl,'location','best');

xlabel('Estimated probability correct')

set(p.h_ax,'ycolor','none','tickdir','out');
xlim([0,1]);

p.next();

POSCLASS = 1;
[X,Y,T,AUC,OPTROCPT,SUBY] = perfcurve(LABELS,SCORES,POSCLASS);

area(X,Y,'EdgeColor','none','FaceColor',0.9*[1,1,1]);
hold all
h = plot(X,Y,'color','k','LineWidth',1);


axis square
hold on
% plot([0,1],[0,1],'k--');
xlabel('P(high confidence | incorrect)');
ylabel('P(high confidence | correct)');


p.format('FontSize',13);
set(h,'LineWidth',2);


end