function run_do_fit_manyseeds()


addpath(genpath('../generic')); 



err_method = 4;

d = load('../data_prepro/LIP_prep_for_DTB/data_for_DTB.mat');

RT  = d.RT;
coh = d.coh;
choice   = d.choice;
group    = d.monkey_id; % monkey id
savefolder = 'fits_neuropixels';


uni_group = unique(group);
nsuj = length(uni_group);
include = ~isnan(RT) & ismember(choice, [0,1]);
ignore_trials = ~include;



%%

dofit_flag = 1;
if dofit_flag

    
%     parpool('local', 2);
    

    % range of the parameters
    kappa = [6,25,17];
    B0    = [0.3,5,1.5];
    a     = [0,4,.1];
    d     = [0,4,1];
    rho   = [-0.7,-0.7,-0.7];
    ndt_mu = [0,0.8,0.35];
    ndt_sigma = [0.01,0.1,0.05];
    coh0 = [-0.1,0.1,0];
    Brectif = [-1,-1,-1];


    params = cat(1,kappa,B0,a,d,rho,ndt_mu,ndt_sigma,coh0,Brectif);
    tl = params(:,1)';
    th = params(:,2)';
    tg = params(:,3)';


    %     Nguess = 10;
    Nguess = 1; % controls if many starting seeds or just 1
    if Nguess==1
        vtg = tg';
    else
        vtg = sample_tguess(tl,th,Nguess,2342344);
    end
    w = cartesian_product(1:nsuj,1:Nguess); % combinations of participant and initial guess

    for i=1:size(w,1)
    % parfor i=1:size(w,1)

        suj_id = w(i,1);
        tg_id  = w(i,2);

        I = group == uni_group(suj_id) & ~ignore_trials; % subset of trials to use

        tg = vtg(:,tg_id)';

        seed = nan;
        pars = struct('suj', suj_id,...
            'sims_per_trial',10,... % how many simulations per trial
            'plotflag',1,... 
            'dt',0.005,... % delta time step
            'err_method', err_method ... % optim criteria
            );
        
        fn_fit = @(theta) (wrapper_DTB_fit(theta,coh(I,:),choice(I),RT(I),seed,pars));

        ptl = tl;
        pth = th;
        [theta,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth);

        tosave = struct('theta',theta,'fval',fval);
        filename = [savefolder,'/fits_',num2str(suj_id),'_',num2str(tg_id)];
        save_parallel(filename,tosave, 0 );
    end


    %% search best
    for i=1:length(uni_group)
        for j=1:Nguess
            filename = [savefolder,'/fits_',num2str(uni_group(i)),'_',num2str(j)];
            aux = load(filename,'theta','fval','tl','tg','th');
            vfval(j) = aux.fval;
        end
        [~,J] = min(vfval);
        filename = [savefolder,'/fits_',num2str(uni_group(i)),'_',num2str(J)];
        aux = load(filename);
        filename_best = [savefolder,'/fits_',num2str(uni_group(i))];
        save_parallel(filename_best, aux, 0 );
    end

end


end
