function H = get_data(dt_ms, neurons_to_include, dataset,varargin)


dt_nans_before_RT = 0.05; % in seconds
positive_is_leftward = 0;
include_mode_config_only = 0;
for i=1:2:length(varargin)
    if isequal(varargin{i},'dt_rel_RT')
        dt_nans_before_RT = varargin{i+1};
    elseif isequal(varargin{i},'positive_is_leftward')
        positive_is_leftward = varargin{i+1};
    elseif isequal(varargin{i},'include_mode_config_only')
        include_mode_config_only = varargin{i+1};
    end
end

extract_neural_data_flag = 1;
if isnan(dt_ms)
    extract_neural_data_flag = 0;
end

if nargin==0 || isempty(dt_ms)
    dt_ms = 5;
end

if nargin<3 || isempty(dataset)
    dataset = 1;
end


if ~isnumeric(dataset)
    datasets = {...
        '20201211_Mars_g0_t2_LIP_4Natalie',...
        '20211011_Jones_2_g0_t1_LIP',...
        '20201030_Mars_g0_t1_LIP',...
        '20201110_Mars_g0_t1_LIP',...
        '20201116_Mars_g0_t1_LIP',...
        '20201208_Mars_g0_t1_LIP',...
        '20211015_Jones_g0_t1_LIP',...
        '20211020_Jones_g0_t0_LIP',...
    };
    dataset = strrep(dataset,'.mat',''); % remove extension
    dataset = find(strcmp(datasets,dataset));
    if isempty(dataset)
        error('invalid dataset');
    end
end

unitIdxLIP_dotsInRF = []; % to update


basedir = fileparts(mfilename('fullpath')); % Get function directory
datadir = fullfile(basedir,'..','data','LIP_Natalie_Gabe2024');


switch dataset
    case 1
        datadir = fullfile(datadir,'20201211');
        filename = '20201211_Mars_g0_t2_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Mars';
        
        
    case 2
        datadir = fullfile(datadir,'20211011');
        filename = '20211011_Jones_2_g0_t1_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Jones';

    case 3
        datadir = fullfile(datadir,'20201030');
        filename = '20201030_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Mars';
        
    case 4
        
        datadir = fullfile(datadir,'20201110');
        filename = '20201110_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Mars';
        
    case 5
        
        datadir = fullfile(datadir,'20201116');
        filename = '20201116_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Mars';
        
    case 6
        
        datadir = fullfile(datadir,'20201208');
        filename = '20201208_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Mars';
        
        
    case 7
        datadir = fullfile(datadir,'20211015');
        filename = '20211015_Jones_g0_t1_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Jones';
        
    case 8
        
        
        datadir = fullfile(datadir,'20211020');
        filename = '20211020_Jones_g0_t0_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        monkey = 'Jones';
        
    
end




%%

if nargin<2 || isempty(neurons_to_include)
    neurons_to_include = 1:length(d(1).spCellPop);
end

if isequal(lower(neurons_to_include),'only_tin')
    % only include Tins?
    neurons_to_include = unitIdxLIP_Tin;
end

%%

trialType = [d.trialType]';

idx_dots = trialType==20; % dots task

dotsOn = [d.dotsOn]';
dotsOff = [d.dotsOff]';
saccadeDetected = [d.saccadeDetected]';

RT = saccadeDetected - dotsOn;
include = idx_dots == 1 & ~isnan(dotsOn) & ~isnan(dotsOff) & ~isnan(RT);

targetsOn = [d.targetsOn]';

coh = [d.sCoh]';
direction = [d.sDir]';
choice = [d.choice]';
correct = [d.correct]';

targ1Pos = cat(1,d.targ1Pos);
targ2Pos = cat(1,d.targ2Pos);


if extract_neural_data_flag
    nTr = length(dotsOn);
    pre_t = 0.2;
    post_t = 5; % 5 seconds max
    nTimeSteps = round(1000*(pre_t+post_t)/dt_ms) + 1;
    nNeurons = length(neurons_to_include);
    H = nan(nNeurons, nTimeSteps, nTr);
    R_ton = nan(nNeurons, nTr);
    for iTr = 1:nTr
        if ~isnan(dotsOn(iTr))
            tini = dotsOn(iTr) - pre_t;
            tend = dotsOn(iTr) + post_t;
            spk = d(iTr).spCellPop(neurons_to_include);
            dt = dt_ms/1000;
            [t,H(:,:,iTr)] = spikeanalysis.spk_to_hist(spk,tini,tend,dt);
            t = t - dotsOn(iTr);
            H(:,(t+dt/2)>[saccadeDetected(iTr) - dotsOn(iTr) - dt_nans_before_RT],iTr) = nan; % NaN spikes after RT

            % target-window
            wt = 0.1; % sec
            [~,aux] = spikeanalysis.spk_to_hist(spk,targetsOn(iTr) - wt, targetsOn(iTr), wt);
            R_ton(:,iTr) = aux(:,1);

        end
    end

    % get the spike times - new
    spike_times_ms = cell(nTr, nNeurons);
    for iTr = 1:nTr
        spk = d(iTr).spCellPop(neurons_to_include);
        for j=1:nNeurons
            spike_times_ms{iTr,j} = 1000*(spk{j}-dotsOn(iTr));
        end
    end



    % extract a higher-res version (1ms), with just the Tin neurons
    dt_high = 1/1000;
    nTimeStepsHigh = [post_t + pre_t]/dt_high + 1;
    DV = nan(nTr, nTimeStepsHigh);
    for iTr = 1:nTr
        if ~isnan(dotsOn(iTr))
            tini = - pre_t;
            tend = + post_t;
            spk = d(iTr).spCellPop(unitIdxLIP_Tin);
            % group them
            spk = {cat(1,spk{:}) - dotsOn(iTr)};
            [td,DV(iTr,:)] = spikeanalysis.spk_to_hist(spk,tini,tend,dt_high);
            %DV(iTr,(td+dt/2)>[saccadeDetected(iTr) - dotsOn(iTr) - 0.05],iTr) = nan; % NaN spikes after RT

        end
    end
    % remove post-decision samples
    DV = motionenergy.remove_post_decision_samples(DV,td, RT); % set nan's after RT
    
    % divide by the number of neurons grouped
    DV = DV / length(unitIdxLIP_Tin);

    DV = DV / (td(2)-td(1)); % for units of firing rate

    %smoothed high
    ns = round(0.08/(td(2)-td(1))); % smoothing time steps for 80ms
    dv_tin_high_smooth = conv2(DV, ones([1,ns]),'same')/ns;
    dv_tin_high_smooth(:,td<-0.1) = nan;
    dv_tin_high_smooth = motionenergy.remove_post_decision_samples(dv_tin_high_smooth, td, RT-dt_nans_before_RT);
    
    

    % get the eye too
    EyeBeforeSacc = [];
    EyeAfterSacc = [];
    for itr = 1:length(d)
        tind1 = d(itr).eye.t>=(d(itr).saccadeDetected-0.1) & d(itr).eye.t<(d(itr).saccadeDetected);
        tind2 = d(itr).eye.t>=(d(itr).saccadeComplete) & d(itr).eye.t<(d(itr).saccadeComplete+0.1);
        
        EyeBeforeSacc = [EyeBeforeSacc; nanmean(d(itr).eye.eyeX(tind1)), nanmean(d(itr).eye.eyeY(tind1))];
        EyeAfterSacc = [EyeAfterSacc; nanmean(d(itr).eye.eyeX(tind2)), nanmean(d(itr).eye.eyeY(tind2))];
        
    end


else
    t = [];
    H = [];
    R_ton = [];
    
end

pulse_on = [d.pulseOn]' - [d.dotsOn]';
pulse_off = [d.pulseOff]' - [d.dotsOn]';
pulse_size = [d.pulseSize]';


coh_extra = coh;
coh_extra(coh==0) = nan;
coh_extra(coh==0 & choice==0) = -0.0000001; % small negative
coh_extra(coh==0 & choice==1) =  0.0000001; % small positive

session = dataset*ones(size(coh));


%%
if (positive_is_leftward) % flip so that positive coherence correspond to a leftward (i.e., contra) choice
    coh = -1*coh;
    coh_extra = -1*coh_extra;
    choice = 1-choice;
    pulse_size = -1*pulse_size;
end
%%
I = include & ~isnan(RT) & [d.complete]'==1;

[~,~,TargConfID] = unique([targ1Pos, targ2Pos],'rows');
isModeTargConfig = TargConfID==mode(TargConfID); % is the target configuration the most common?

if include_mode_config_only==1
    I = I & isModeTargConfig==1;
end


monkey = repmat({monkey}, length(coh), 1);
if extract_neural_data_flag
    out = struct('H',H(:,:,I), 'unitIdxLIP_Tin',unitIdxLIP_Tin, 'unitIdxLIP_Tout', unitIdxLIP_Tout, ...
        't',t, 'RT',RT(I), 'included',I, 'coh',coh(I), 'coh_extra', coh_extra(I), ...
        'direction',direction(I), 'choice',choice(I), ...
        'correct',correct(I), 'unitIdxLIP_dotsInRF',unitIdxLIP_dotsInRF, 'pulse_on',pulse_on(I), ...
        'targ1Pos',targ1Pos(I,:),'targ2Pos',targ2Pos(I,:),'isModeTargConfig',isModeTargConfig(I),...
        'pulse_off',pulse_off(I), 'pulse_size',pulse_size(I),'dataset',filename,'idataset',dataset,...
        'R_ton',R_ton(:,I),'td',td,'DV',DV(I,:),'DV_smooth',dv_tin_high_smooth(I,:),...
        'nneurons',size(H,1), 'ntimes',size(H,2),'ntrials',sum(I),'session',session(I),...
        'eyeBeforeSacc',EyeBeforeSacc(I,:),'eyeAfterSacc',EyeAfterSacc(I,:));

    out.monkey = monkey(I);
    out.spike_times_ms = spike_times_ms(I,:);

else
    out = struct('RT',RT(I), 'included',I, 'coh',coh(I), 'coh_extra',coh_extra(I), 'drection',direction(I), 'choice',choice(I), ...
        'correct',correct(I), 'pulse_on',pulse_on(I), ...
        'pulse_off',pulse_off(I), 'pulse_size',pulse_size(I),'dataset',filename,'idataset',dataset);

    out.monkey = monkey(I);
end


file_new_tin = fullfile(basedir,'prepro_recalc_Tin_Tout',['d',num2str(dataset),'.mat']);
if extract_neural_data_flag && exist(file_new_tin,'file')
    aux = load(file_new_tin);
    out.idx_Tin = aux.idx_Tin;
    out.idx_Tout = aux.idx_Tout;
end

% if extract_neural_data_flag
%     out.minCells = minCells;
% end

H = out; % first output




end

