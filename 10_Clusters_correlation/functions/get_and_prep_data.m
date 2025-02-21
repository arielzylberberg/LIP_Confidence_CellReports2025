function S = get_and_prep_data(win)

if nargin==0 || isempty(win)
    win = 0.05;
end

load('../03_Kmeans_calculations/output_data/data_per_cluster.mat','S');


%% 
for i=1:length(S)


    %% convolve
    dt = S(i).t(2) - S(i).t(1);
    sm = round(win/dt); 
    H = conv2(1,ones(sm,1),S(i).H,'same')/sm; % H has the spike counts grouping neurons belonging to the same cluster

    %% remove last 50 ms
    H = motionenergy.remove_post_decision_samples(H, S(i).t, S(i).RT-0.05); 

    %% standarize`
    zHres = nanzscore_bygroup(H, [S(i).coh, S(i).dataset]);

    %% remove baseline from each trial
    baseline_subtract = 0;
    if baseline_subtract
        tind = S(i).t>=-0.1 & S(i).t<=0;
        zHres = zHres - nanmean(zHres(:,tind),2);
    end

    %%
    S(i).zHres = zHres;

    %% sub-sample
    step = round(sm);
    S(i).t = S(i).t(1:step:end);
    S(i).H = S(i).H(:,1:step:end);
    S(i).zHres = S(i).zHres(:,1:step:end);

end
