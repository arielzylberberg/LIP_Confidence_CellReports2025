function out = get_data_passive_motion(dt_ms, neurons_to_include, dataset,varargin)


if nargin==0 || isempty(dt_ms)
    dt_ms = 5;
end

if nargin<3 || isempty(dataset)
    dataset = 1;
end


positive_is_leftward = 0;
for i=1:2:length(varargin)
    if isequal(varargin{i},'positive_is_leftward')
        positive_is_leftward = varargin{i+1};
    end
end

% datasets = {'20201211_Mars','Jones_20211011_4Natalie','20201030_Mars_g0_t1_LIP','20201110_Mars_g0_t1_LIP','20201116_Mars_g0_t1_LIP','20201208_Mars_g0_t1_LIP'};

if ~isnumeric(dataset)
    datasets = {'20201211_Mars_g0_t2_LIP_4Natalie','20211011_Jones_2_g0_t1_LIP','20201030_Mars_g0_t1_LIP',...
        '20201110_Mars_g0_t1_LIP','20201116_Mars_g0_t1_LIP','20201208_Mars_g0_t1_LIP',...
        '20211015_Jones_g0_t1_LIP','20211020_Jones_g0_t0_LIP','20211007_Neo','20211004_Neo'};
    dataset = strrep(dataset,'.mat',''); % remove extension
    dataset = find(strcmp(datasets,dataset));
    if isempty(dataset)
        error('invalid dataset');
    end
end




% basedir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/93 - LIP_neuropixel_dots';
basedir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels';
datadir = fullfile(basedir,'data');

% load('../data/minCells.mat');
aux = load(fullfile(datadir,'newNeuronClasses.mat'));
cells = aux.newNeuronClasses;
clear cell_dates
for i=1:length(cells)
    cell_dates{i} = cells{i}.date;
end

switch dataset
    case 1
        datadir = fullfile(datadir,'20201211_Mars');
        filename = '20201211_Mars_g0_t2_LIP_4Natalie';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         unitIdxLIP_dotsInRF = [];
        
        I = ismember(cell_dates,'201211M');
        cells = cells{I};
        
        
    case 2
        datadir = fullfile(datadir,'Jones_20211011_4Natalie');
        filename = '20211011_Jones_2_g0_t1_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         load(fullfile(datadir,'unitIdxLIP_dotsInRF.mat'));
        
        I = ismember(cell_dates,'211011J');
        cells = cells{I};

    case 3
        datadir = fullfile(datadir,'20201030_Mars_g0_t1_LIP');
        filename = '20201030_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         unitIdxLIP_dotsInRF = [];
        
        I = ismember(cell_dates,'201030M');
        cells = cells{I};
        
    case 4
        
        datadir = fullfile(datadir,'20201110_Mars_g0_t1_LIP');
        filename = '20201110_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         unitIdxLIP_dotsInRF = [];
        
        I = ismember(cell_dates,'201110M');
        cells = cells{I};
        
    case 5
        
        datadir = fullfile(datadir,'20201116_Mars_g0_t1_LIP');
        filename = '20201116_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         unitIdxLIP_dotsInRF = [];
        
        I = ismember(cell_dates,'201116M');
        cells = cells{I};
        
    case 6
        
        datadir = fullfile(datadir,'20201208_Mars_g0_t1_LIP');
        filename = '20201208_Mars_g0_t1_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         unitIdxLIP_dotsInRF = [];
        
        I = ismember(cell_dates,'201208M');
        cells = cells{I};
        
        
    case 7
        datadir = fullfile(datadir,'20211015_Jones');
        filename = '20211015_Jones_g0_t1_LIP';
        load(fullfile(datadir,filename));
%         unitIdxLIP_Tout = [];
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         load(fullfile(datadir,'unitIdxLIP_dotsInRF.mat'));
        
        I = ismember(cell_dates,'211015J');
        cells = cells{I};
        
    case 8
        
        
        datadir = fullfile(datadir,'20211020_Jones');
        filename = '20211020_Jones_g0_t0_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         load(fullfile(datadir,'unitIdxLIP_dotsInRF.mat'));
        
        I = ismember(cell_dates,'211020J');
        cells = cells{I};
        
    case 9
        
        datadir = fullfile(datadir,'20211007_Neo');
        filename = 'Neo211007_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
%         aux = load(fullfile(datadir,'unitIdxLIP_dotsRF.mat'));
%         unitIdxLIP_dotsInRF = aux.unitIdxLIP_dotsRF;
        
        I = ismember(cell_dates,'211007N');
        if sum(I)>0
            cells = cells{I};
        else
            cells = [];
        end
        
    case 10
        
        datadir = fullfile(datadir,'20211004_Neo');
        filename = 'Neo211004_LIP';
        load(fullfile(datadir,filename));
%         load(fullfile(datadir,'unitIdxLIP_Tout.mat'));
%         load(fullfile(datadir,'unitIdxLIP_Tin.mat'));
        unitIdxLIP_dotsInRF = [];
        
        I = ismember(cell_dates,'211004N');
        if sum(I)>0
            cells = cells{I};
        else
            cells = [];
        end
        
end

unitIdxLIP_Tin = cells.unitIdxLIP_Tin;
unitIdxLIP_Tout  = cells.unitIdxLIP_Tout;
unitIdxLIP_dotsInRF = [];
minCells = struct('DinRFc', cells.unitIdx_DinRFcC ,'DinRFi',cells.unitIdx_DinRFiC);



%%

if nargin<2 || isempty(neurons_to_include)
    neurons_to_include = 1:length(d(1).spCellPop);
end

%%

trialType = [d.trialType]';

idx_viewing = trialType==50; % dots-viewing task

dotsOn = [d.dotsOn]';
dotsOff = [d.dotsOff]';

include = idx_viewing == 1 & ~isnan(dotsOn) & ~isnan(dotsOff);


coh = [d.sCoh]';
direction = [d.sDir]';
choice = [d.choice]';
correct = [d.correct]';


nTr = length(dotsOn);
pre_t = 0.2;
post_t = 3; % 3 seconds max
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
        H(:,(t+dt/2)>[dotsOff(iTr) - dotsOn(iTr)],iTr) = nan; % NaN spikes after dots-off

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

%%
if (positive_is_leftward)
    choice = 1-choice;
    coh = -1*coh;
end

%%

I = include & [d.complete]'==1;

out = struct('H',H(:,:,I), 'unitIdxLIP_Tin',unitIdxLIP_Tin, 'unitIdxLIP_Tout', unitIdxLIP_Tout, ...
    't',t, 'included',I, 'coh',coh(I), 'direction',direction(I), ...
    'unitIdxLIP_dotsInRF',unitIdxLIP_dotsInRF, ...
    'dataset',filename,...
    'dot_dur',dotsOff(I)-dotsOn(I),...
    'nneurons',size(H,1), 'ntimes',size(H,2),'ntrials',sum(I));

out.spike_times_ms = spike_times_ms(I,:);


file_new_tin = ['../prepro_recalc_Tin_Tout/d',num2str(dataset),'.mat'];
if exist(file_new_tin,'file')
    aux = load(file_new_tin);
    out.idx_Tin = aux.idx_Tin;
    out.idx_Tout = aux.idx_Tout;
end

out.minCells = minCells;



end

