function [Bopt,RTopt,AccOpt] = optim_bound_single_coh(coh,kappa,tnd,tinter,sigma)
% [Bopt,RTopt,AccOpt] = optim_bound_single_coh(coh,kappa,tnd,tinter,sigma)
% example: [Bopt,RTopt,AccOpt] = dotsanalysis.optim_bound_single_coh(.1,12,.3,3,sqrt(1.2))

%% optimal av rew

Pplus  = inline('1./(1+exp(-(2*k*C*A/sigma^2)))','k','C','A','sigma');
RTmean = inline('A./(k*C).*tanh(k*C*A/sigma^2)','k','C','A','sigma');

%
% sigma = 1;
C = coh;
k = kappa;
vA = linspace(0,3,1000);
iti = tnd + tinter; %inter trial interval
% figure
for i=1:length(vA)
    acc = Pplus(k,C,vA(i),sigma);
    rt(i) = RTmean(k,C,vA(i),sigma);
    av_rew(i) = acc / (rt(i)+iti);
end

% plot(vA,av_rew)
Bopt = vA(av_rew==max(av_rew));
RTopt = rt(av_rew==max(av_rew));
AccOpt = Pplus(k,C,Bopt,sigma);


% % 
% %% compare with DTB, for different dt
% coh = [-C,C]';
% rt = nan(size(coh));
% choice = nan(size(coh));
% c = nan(size(coh));
% plot_flag=0;
% kappa = k;
% theta = [kappa,nan,nan,Bopt,0,0,0,0];
% 
% vDT = [0.002,0.001,0.0005,0.0001,0.00001];
% for i=1:length(vDT)
%     pars.methodforward_flag = 1;
%     dt = vDT(i);
%     pars.t = [0:dt:10];
%     [err,P] = wrapper_dtb_parametricbound_rt_2(theta,rt,coh,choice,c,pars,plot_flag);
%     RT_DTB(i) = P.up.mean_t(1);
%     Acc_DTB(i) = P.up.p(2);
% end
% 
% %%
% p = publish_plot(1,2);
% p.next();
% plot(1000*vDT,RT_DTB,'o-')
% hold all
% plot(xlim,[RTopt,RTopt],'r')
% xlabel('Time step in DTB (ms)')
% ylabel('RT')
% 
% p.next();
% plot(1000*vDT,Acc_DTB,'o-')
% hold all
% plot(xlim,[AccOpt,AccOpt],'r')
% xlabel('Time step in DTB (ms)')
% ylabel('Acc')
% 
% hl = legend('DTB','theory');
% set(hl,'location','best')
% p.format()
% set(gcf,'Position',[369  153  821  374])
% 
% export_fig('-eps','fig_DTB_accuracy_vs_time_step')