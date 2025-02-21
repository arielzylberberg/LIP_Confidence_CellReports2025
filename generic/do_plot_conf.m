function p = do_plot_conf(conf, coh, RT, correct, varargin)

npoints = 300;
nbins = 5;
for i=1:length(varargin)
    if isequal(varargin{i},'npoints')
        npoints = varargin{i+1};
    end
end



plotflag = 1;
filt = [];

duration = RT;
scoh = abs(coh);

p = publish_plot(1,3);
set(gcf,'Position',[298  415  802  215]);
p.next();
min_tr_per_cond = 4;
[tt,xx,ss] = curva_media(conf, 100*scoh, correct==1,0);
[~,ntr_per_cond] = curva_suma(ones(size(conf)), 100*scoh, correct==1,0);
ind = ntr_per_cond>=min_tr_per_cond;
terrorbar(tt(ind),xx(ind),ss(ind),'color','b','marker','o','markerfacecolor','b');
hold all
[tt,xx,ss] = curva_media(conf, 100*scoh, correct==0,0);
[~,ntr_per_cond] = curva_suma(ones(size(conf)), 100*scoh, correct==0,0);
ind = ntr_per_cond>=min_tr_per_cond;
terrorbar(tt(ind),xx(ind),ss(ind),'color','r','marker','o','markerfacecolor','r');
% hl = legend('correct','incorrect');
xlabel('Motion strength [%]');
ylabel('Confidence estimate');
% set(hl,'location','southeast');


p.next();
[X,Y] = dotsanalysis.choice_against_duration_moving(conf,correct==0,duration,npoints,0,filt);

plot(X{1},Y{1},'b');
hold all;
plot(X{2},Y{2},'r');
% [X,Y,E] = dotsanalysis.choice_against_duration_var(d,scoh,duration,nbins,plotflag,filt);
xlabel('Reaction time [s]');
% ylabel('Confidence estimate');
% hl = legend('correct','error');

% p.format('FontSize',18,'LineWidthPlot',1);

% set(hl,'location','best');


p.next();
[X,Y] = dotsanalysis.choice_against_duration_moving(conf,scoh,duration,npoints, 0 ,correct==1);
color = movshon_colors(6);
for i=1:length(X)
    plot(X{i},Y{i},'color',color(i,:));
    hold all
end
% [X,Y,E] = dotsanalysis.choice_against_duration_var(conf,scoh,duration,nbins,plotflag,filt);

xlabel('Reaction time [s]');
% ylabel('Decoders conf in choice');

% hl = legend_n(unique(100*scoh));


% ylabel('Decoders conf in choice');

p.format('FontSize',12,'LineWidthPlot',1.0,'MarkerSize',6);
% set(hl,'FontSize',11);
% set(hl,'Position',[0.8830 0.5411 0.0779 0.3930]);

% same_ylim(p.h_ax);
% set(p.h_ax,'ylim',[0.6,1]);

set(p.h_ax,'ylim',[0.5,1]);

end
