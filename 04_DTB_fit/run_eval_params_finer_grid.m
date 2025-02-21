addpath(genpath('../generic'));

%%

d = load('../data_prepro/LIP_prep_for_DTB/data_for_DTB.mat');
RT = d.RT;
coh = d.coh;
choice = d.choice;
group = d.monkey_id;
datafolder = 'fits_neuropixels';


uni_group = unique(group);
nsuj = length(uni_group);


%%

colores = [173,93,0; 0,0,0]/255;

p = publish_plot(2,1);
set(gcf,'Position',[663  488  287  490]);

for suj = 1:2
    include = ~isnan(RT) & ismember(choice, [0,1]) & group==suj;
    ignore_trials = ~include;
    
    aux = load([datafolder,'/fits_',num2str(suj),'.mat']);
    theta = aux.theta;
    
    I = include;
    
    seed = 447478 + suj;
    uni_coh_fine = linspace(-0.512,0.512,51)';
    pars = struct('sims_per_trial', 4000, 'plotflag',0,'suj',unique(group(I)),'dt',0.005);
    NaNs = nan(size(uni_coh_fine));
    [~,m] =  wrapper_DTB_fit(theta,uni_coh_fine,NaNs,NaNs,seed,pars);


    %% plot choice and RTs
    
    p.next();
    [tt,xx,ss] = curva_media(RT, 100*coh,include,0);
    terrorbar(tt,xx,ss,'color',colores(suj,:),'marker','.','LineStyle','none','markerfacecolor',colores(suj,:),'markeredgecolor',...
        colores(suj,:));
    hold all
    [tt,xx] = curva_media(m.RT, 100*m.coh,[],0);
    h(suj) = plot(tt,xx,'color',colores(suj,:));
    
    ylabel('Response time [s]');
%     xlabel('Motion strength')
    p.format();

    p.next();
    m_choice = m.winner==2; % left
    [tt,xx,ss] = curva_media(choice, 100*coh,include,0);
    terrorbar(tt,xx,ss,'color',colores(suj,:),'marker','.','LineStyle','none','markerfacecolor',colores(suj,:),'markeredgecolor',...
        colores(suj,:));
    hold all
    [tt,xx] = curva_media(m_choice, 100*m.coh,[],0);
    plot(tt,xx,'color',colores(suj,:));
    ylabel('Prop. leftward choice');
    xlabel('Motion strength [%coh]');


end

p.current_ax(1);
hl = legend(h,'Monkey M','Monkey J');
set(hl,'location','best','box','off');

% p.format('MarkerSize',[20,20]);
p.format('MarkerSize',[17,17],'LineWidthPlot',0.75);

p.saveas('./figures/fig_fits_with_lower_bound');


% p.append_to_pdf('fig_psych',1,1);


