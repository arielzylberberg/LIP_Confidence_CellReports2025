% load ../analysis_kmeans/out_3clusters.mat
addpath('../generic'); 

%%
ndatasets = 8;

dataset = [];
choice = [];
RT = [];
correct = [];
coh = [];
coh_extra = [];
H = [];
neuron_id = [];
spk_RT = [];
spk_suma = []; 
trnum = [];

for idataset = 1:ndatasets

    disp(num2str(idataset));

    dt_ms = 200;
    d = get_data(dt_ms,'only_Tin',idataset,'positive_is_leftward',1);

    nneurons = size(d.spike_times_ms,2);

    spk_RT_aux = count_spikes_in_window(d.spike_times_ms, [], 1000*(d.RT-0.15), 1000*(d.RT-0.05));
    spk_RT = [spk_RT; spk_RT_aux(:)];


    spk_suma_aux = count_spikes_in_window(d.spike_times_ms, [], 1000*( 0.2 ), 1000*(d.RT-0.05));
    spk_suma = [spk_suma; spk_suma_aux(:)];


    for j = 1:nneurons

        deltat_ms = 5;
        [t,Haux] = spikeanalysis.spk_to_hist(d.spike_times_ms(:,j) , -200, 4000, deltat_ms);
    

        ntr = length(d.choice);
    
        dataset = [dataset; idataset(ones(ntr,1))];
        choice  = [choice; d.choice];
        RT      = [RT; d.RT];
        correct = [correct; d.correct];
        coh     = [coh; d.coh];
        coh_extra = [coh_extra; d.coh_extra];
        H       = [H; Haux];
%         cluster = [cluster; ones(ntr,1)*clust(j)];
        neuron_id = [neuron_id; ones(ntr,1)*j];
        trnum = [trnum; [1:ntr]'];
    end

end


[~,~,neuron_id] = unique([dataset, neuron_id],'rows');

t = t/1000;
save('data_Tin','dataset','choice','RT','correct','coh','coh_extra','spk_RT','spk_suma','neuron_id','trnum'); % no H
save('data_Tin_with_H','dataset','choice','RT','correct','coh','coh_extra','spk_RT','spk_suma','neuron_id','t','trnum','H','-v7.3');

