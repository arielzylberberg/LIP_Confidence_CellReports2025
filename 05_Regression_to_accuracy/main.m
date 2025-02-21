addpath('../generic'); 

%%
rng(19858,'twister');

labels = {
    'TinRateRT',
    'TinMeanRT',
    'TinCluster1RT',
    'TinCluster2RT',
    'TinCluster3RT',
    'TinCluster13RT',
    'TinCluster13AverageRT',
    'EyeBeforeAfter', % for reviewers
    'TinRateRT_shorter_window' % for reviewers
};


regression_to_conf_group_sessions(labels);