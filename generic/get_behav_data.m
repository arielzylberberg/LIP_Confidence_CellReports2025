function [RT, coh, direction, choice, correct, data_set, ...
    pulse_on, pulse_off, pulse_size] = get_behav_data()


coh = [];
direction = [];
choice = [];
correct = [];
pulse_on = [];
pulse_off = [];
pulse_size = [];
data_set = [];
RT = [];


v_dataset = 1:6; % all

for idata = 1:length(v_dataset)
    dataset = v_dataset(idata);
    
    switch dataset
        case 1
            datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201211_Mars';
            load(fullfile(datadir,'20201211_Mars_g0_t2_LIP_4Natalie'));
            
            
        case 2
            datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/Jones_20211011_4Natalie';
            load(fullfile(datadir,'20211011_Jones_2_g0_t1_LIP'));
            
            
        case 3
            datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201030_Mars_g0_t1_LIP';
            load(fullfile(datadir,'20201030_Mars_g0_t1_LIP'));
            
        case 4
            
            datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201110_Mars_g0_t1_LIP';
            load(fullfile(datadir,'20201110_Mars_g0_t1_LIP'));
            
        case 5
            
            datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201116_Mars_g0_t1_LIP';
            load(fullfile(datadir,'20201116_Mars_g0_t1_LIP'));
            
        case 6
            
            datadir = '/Users/arielzy/Dropbox/Data/1 - LAB/01 - Proyectos/96 - LIP_Confidence_neuropixels/data/20201208_Mars_g0_t1_LIP';
            load(fullfile(datadir,'20201208_Mars_g0_t1_LIP'));
            
            
    end
    
    
    
    
    %%
    
    trialType = [d.trialType]';
    
    idx_dots = trialType==20; % dots task
    
    dotsOn = [d.dotsOn]';
    dotsOff = [d.dotsOff]';
    saccadeDetected = [d.saccadeDetected]';
    
    RTx = saccadeDetected - dotsOn;
    
    include = idx_dots == 1 & ~isnan(dotsOn) & ~isnan(dotsOff) & ~isnan(RTx);
    RT = [RT; RTx(include)];
    
    cohx = [d.sCoh]';
    directionx = [d.sDir]';
    choicex = [d.choice]';
    correctx = [d.correct]';
    
    
    coh = [coh; cohx(include)];
    direction = [direction; directionx(include)];
    choice = [choice; choicex(include)];
    correct = [correct; correctx(include)];
    data_set = [data_set; ones(sum(include),1)*dataset];
    
    
    pulse_on = [pulse_on; [d(include).pulseOn]' - [d(include).dotsOn]'];
    pulse_off = [pulse_off; [d(include).pulseOff]' - [d(include).dotsOn]'];
    pulse_size = [pulse_size; [d(include).pulseSize]'];
    
    
end

