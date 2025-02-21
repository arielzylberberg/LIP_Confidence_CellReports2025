function do_save_data_for_kmeans_analysis()


addpath('../generic/');


%%

C = load('./output_data/out_3clusters.mat');


ndatasets = 8;


str = {'all Tin','Cluster 1','Cluster 2','Cluster 3'};
nclusters = 3;

for i=1:length(str)
    S(i).dataset = [];
    S(i).choice = [];
    S(i).RT = [];
    S(i).correct = [];
    S(i).coh = [];
    S(i).H = [];
    S(i).pulse_on = [];
    S(i).pulse_off = [];
    S(i).pulse_size = [];
    S(i).str = str{i};
    S(i).spk_RT = [];
end


for idataset = 1:ndatasets

    disp(num2str(idataset));

    dt_ms = 200;
    d = get_data(dt_ms,[],idataset,'positive_is_leftward',1);
    
    clust = C.clusters(C.data_id==idataset);
    ntr = length(d.RT);
    for i=1:length(str)
        switch i
            case 1
                I = ismember(clust, 1:nclusters);
                I = d.unitIdxLIP_Tin(I);
            case 2
                I = ismember(clust, 1);
                I = d.unitIdxLIP_Tin(I);
            case 3
                I = ismember(clust, 2);
                I = d.unitIdxLIP_Tin(I);
            case 4
                I = ismember(clust, 3);
                I = d.unitIdxLIP_Tin(I);

        end
        s = d.spike_times_ms(:,I);
        spk = cell(ntr,1);
        for j=1:ntr
            spk{j} = sort(cat(1,s{j,:}));
        end

        deltat_ms = 1;
        
        [t,Haux] = spikeanalysis.spk_to_hist(spk , -400, 5000, deltat_ms);
        
        nneurons = length(I);
        Haux = Haux/nneurons; % to make it per neuron


        S(i).dataset    = [S(i).dataset; idataset(ones(ntr,1))];
        S(i).choice     = [S(i).choice; d.choice];
        S(i).RT         = [S(i).RT; d.RT];
        S(i).correct    = [S(i).correct; d.correct];
        S(i).coh        = [S(i).coh; d.coh];
        S(i).H          = [S(i).H; Haux];
        S(i).pulse_on   = [S(i).pulse_on; d.pulse_on];
        S(i).pulse_off  = [S(i).pulse_off; d.pulse_off];
        S(i).pulse_size     = [S(i).pulse_size; d.pulse_size];
        

        S(i).t = t/1000;
        
        if numel(I)==0
            aux = nan(size(d.choice));
        else
            aux = count_spikes_in_window(d.spike_times_ms, I, d.RT*1000-150, d.RT*1000-50)';
            if size(aux,1)>1
                aux = nanmean(aux)';
            else
                aux = aux';
            end
        end
        S(i).spk_RT = [S(i).spk_RT; aux];

    end
    
end



t = t/1000;
save('./output_data/data_per_cluster','S','-v7.3'); 
