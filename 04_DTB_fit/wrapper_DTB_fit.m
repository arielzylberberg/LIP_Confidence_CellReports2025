function [err,m] = wrapper_DTB_fit(theta,coh,choices,RT,seed,pars)

kappa   = theta(1);
B0      = theta(2);
a       = theta(3);
d       = theta(4);
rho     = theta(5);
ndt_mu  = theta(6);
ndt_sigma = theta(7);
coh0    = theta(8);

if length(theta)>8
    Brectif = theta(9);
else
    Brectif = 0;
end


%%

sims_per_trial = pars.sims_per_trial;


ntr = length(coh);
ind = 1:ntr;
tr_id = repmat(ind,sims_per_trial,1);
tr_id = tr_id(:);


delta_t = 0.01;
if isfield(pars,'dt')
    delta_t = pars.dt;
end



tau_kernel_noise = 0;
tdelay_acum = 0.0; 
y0std = 0;

% this does the simulations...
m = dtb_fake_data('Brectif',Brectif,'srho', rho ,'USfunc','Linear','a',a,...
    'd',d + tdelay_acum,'t',[0:delta_t:10],'B0',B0,'kappa',kappa,...
    'y0',0 ,...
    'y0std',y0std ,... % testing
    'coh0',coh0,...
    'coh',coh(tr_id,:),...
    'ndt_mu',ndt_mu, ...
    'ndt_sigma',ndt_sigma,...
    'ndt_method_flag', 2,...
    'boost',1,...
    'attention_shift_rate', 0,...
    'tau_kernel_noise', tau_kernel_noise,... % seconds- kernel to smooth the noise
    'tdelay_acum', tdelay_acum,... % seconds - delay between the ev stream and the acum
    'tdelay_acum_std', 0,...
    'urgency_flag', 1, ... % models bound as urgency
    'seed',seed);

m.make_fakedata();
m.diffuse_to_bound();


%% compute error according to some criteria (err_method)

dt = m.t(2)-m.t(1);
        

% err method: max likelihood
p_RT_choice = nan(size(RT,1),1);
uni_coh = nanunique(coh);
n = length(uni_coh);
filt = RT>ndt_mu;
for i=1:n % unique coh
    for k=1:2 % choice
        K = coh==uni_coh(i) & filt;
        I = K & choices==(k-1);
        J = m.coh == uni_coh(i) & m.winner==k;
        if sum(J)>1 && sum(I)>0 && sum(~isnan(m.decision_time(J)))>1
            pd = fitdist(m.decision_time(J),'kernel','Kernel','epanechnikov','support','positive');

            pdf = pd.pdf(m.t);
            pdf = pdf/sum(pdf);
            pd_ndt = makedist('Normal','mu',ndt_mu,'sigma',ndt_sigma); % non-dec time distribution
            pd_trunc = truncate(pd_ndt,0,inf); % truncate it
            ndt = pd_trunc.pdf(m.t)*dt;
            pdf = conv(ndt, pdf); % convolve with the decision time distribution
            pdf = pdf(1:length(m.t));
            dt = m.t(2)-m.t(1);

            rt_step = ceil(RT(I)/dt);
            p_RT_choice(I) = pdf(rt_step) * nanmean(choices(K)==(k-1)); % p(RT|choice) * p(choice) = p (RT,choice)
        end
    end
end

pPred = p_RT_choice;
pPred(pPred<=0 | isnan(pPred)) = eps;
err = -sum(log(pPred));


if pars.plotflag
    figure(pars.suj)
    clf
    subplot(1,2,1)
    curva_media(choices,coh,[],1);
    hold on
    curva_media(m.winner-1,m.coh,[],1);
    xlabel('Motion Coherence');
    ylabel('Choice');

    subplot(1,2,2)
    curva_media(RT,coh,[],1);
    hold on
    curva_media(m.RT,m.coh,[],1);
    ylabel('RT [s]');
    xlabel('Choice');

    drawnow
    
end


%% print
var_names = {'k','B0','a','d','rho','ndt_m','ndt_s','coh0','Brectif'};
fprintf_params(var_names, err, [kappa, B0, a, d, rho, ndt_mu, ndt_sigma, coh0, Brectif]);


end

