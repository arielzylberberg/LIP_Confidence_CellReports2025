function [H, unitIdxLIP_Tin, unitIdxLIP_Tout, t, RT, L, target_pos, ...
    uni_target_pos, mL, target_pos_dots_task, AUC] = ...
    get_data_memory_saccades(dt_ms, dataset)


if nargin==0 || isempty(dt_ms)
    dt_ms = 5;
end

if isempty(dataset) || nargin<2
    dataset = 1;
end

% datasets = {'20201211_Mars','Jones_20211011_4Natalie','20201030_Mars_g0_t1_LIP','20201110_Mars_g0_t1_LIP','20201116_Mars_g0_t1_LIP','20201208_Mars_g0_t1_LIP'};

switch dataset
    case 1
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201211_Mars';
        load(fullfile(datadir,'20201211_Mars_g0_t2_LIP_4Natalie'));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
    case 2
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/Jones_20211011_4Natalie';
        load(fullfile(datadir,'20211011_Jones_2_g0_t1_LIP'));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        load(fullfile(datadir,'unitIdxLIP_dotsInRF.mat'));

    case 3
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201030_Mars_g0_t1_LIP';
        load(fullfile(datadir,'20201030_Mars_g0_t1_LIP'));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
    case 4
        
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201110_Mars_g0_t1_LIP';
        load(fullfile(datadir,'20201110_Mars_g0_t1_LIP'));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
    case 5
        
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201116_Mars_g0_t1_LIP';
        load(fullfile(datadir,'20201116_Mars_g0_t1_LIP'));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
    case 6
        
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201208_Mars_g0_t1_LIP';
        load(fullfile(datadir,'20201208_Mars_g0_t1_LIP'));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
    case 7
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20211015_Jones';
        filename = '20211015_Jones_g0_t1_LIP';
        load(fullfile(datadir,filename));
        unitIdxLIP_Tout = [];
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        load(fullfile(datadir,'unitIdxLIP_dotsInRF.mat'));
        
        
    case 8
        
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20211020_Jones';
        filename = '20211020_Jones_g0_t0_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        load(fullfile(datadir,'unitIdxLIP_dotsInRF.mat'));
        
        
    case 9
        
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20211007_Neo';
        filename = 'Neo211007_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        aux = load(fullfile(datadir,'unitIdxLIP_dotsRF.mat'));
        unitIdxLIP_dotsInRF = aux.unitIdxLIP_dotsRF;
        
        
    case 10
        
        datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20211004_Neo';
        filename = 'Neo211004_LIP';
        load(fullfile(datadir,filename));
        load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
        load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
        
        
end


%%



%%

trialType = [d.trialType]';

idx_mem_sacc = trialType==3; % memory saccades

targetsOn = [d.targetsOn]';
fixOff = [d.fixOff]';
saccadeDetected = [d.saccadeDetected]';

RT = fixOff - targetsOn; % not really RT
correct = [d.correct]';
t_ignore_init = 0.2;
include = idx_mem_sacc == 1 & ~isnan(targetsOn)  & ~isnan(RT) & [d.complete]'==1 & ...
    correct==1 & RT>t_ignore_init;


target_pos = cat(1,d.targ1Pos);


nTr = length(targetsOn);
base_t = -0.2;
max_t = 3;
nTimeSteps = round(1000*(max_t - base_t)/dt_ms) + 1;
nNeurons = length(d(1).spCellPop);
H = nan(nNeurons, nTimeSteps, nTr);
L = nan(nNeurons, nTr);
for iTr = 1:nTr
    if include(iTr)==1
        tini = targetsOn(iTr) + base_t;
        tend = targetsOn(iTr) + max_t;
        spk = d(iTr).spCellPop;
        dt = dt_ms/1000;
        [t,H(:,:,iTr)] = spikeanalysis.spk_to_hist(spk,tini,tend,dt);
        t = t - targetsOn(iTr);
        % H(:,(t+dt/2)>[saccadeDetected(iTr) - targetsOn(iTr)],iTr) = nan; % NaN spikes after RT
        H(:,(t+dt/2)>[fixOff(iTr) - targetsOn(iTr)],iTr) = nan; % NaN spikes after RT
        
        % ignore the first 200ms and count spikes
        
        tind = t>t_ignore_init;
        L(:, iTr) = squeeze(nansum(H(:, tind,iTr),2));
    end
end

%%
H = H(:,:,include);
RT = RT(include);
target_pos = target_pos(include,:);
L = L(:,include);

L = bsxfun(@times, L', 1./(RT - t_ignore_init)); % convert to spikes per sec

% now average L per unique value of target_pos
[uni_target_pos,~,idx] = unique(target_pos,'rows');
[~,mL] = curva_media(L,idx,[],0);
mL = mL';


if 0 % as it was until May 21st, 2024
    I = find(trialType == 20,1); % dots
    
    target_pos_dots_task.targ1 = d(I).targ1Pos;
    target_pos_dots_task.targ2 = d(I).targ2Pos;

else
    % alternative: most freq. config. 
    I = trialType==20; % dots
    [uni_targ_pos,~,count] = unique([cat(1,d(I).targ1Pos), cat(1,d(I).targ2Pos)],'rows');
    target_pos_dots_task.targ1 = uni_targ_pos(mode(count),1:2);
    target_pos_dots_task.targ2 = uni_targ_pos(mode(count),3:4);

end
%% calc dprime for the two relevant targets
T = zeros(size(H,3),1);

I = all(target_pos==target_pos_dots_task.targ1,2);
T(I) = 1;
I = all(target_pos==target_pos_dots_task.targ2,2);
T(I) = 2;

AUC = nan(nNeurons,1);
if sum(T==1)>0 && sum(T==2)>0
    for ineuron = 1:nNeurons
        labels = [ones(sum(T==1),1); 2*ones(sum(T==2),1)];
        scores = [L(T==1,ineuron);L(T==2,ineuron)];
    
        [X,Y,~,AUC(ineuron)] = perfcurve(labels,scores,2);
    end
end



end

