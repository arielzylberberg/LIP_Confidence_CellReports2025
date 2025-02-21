


addpath('../generic/');

%%

num_datasets = 8;
cont = 0;
clear B

for idataset=1:num_datasets
    d = get_data(1,'only_tin',idataset,'dt_rel_RT',0);
    for j=1:d.nneurons
        h = squeeze(d.H(j,:,:));
    

        H = motionenergy.remove_post_decision_samples(h', d.t, d.RT-0.05);
        [datamp,timeslocked] = eventlockedmatc(H, d.t, d.pulse_on);
    
        pulse_size = d.pulse_size;
        coh = d.coh;
    
        ts = [-60:20:800]/1000;
        cont = cont+1;
        for i=1:length(ts)
    
            tind = timeslocked>=ts(i) & timeslocked<=(ts(i) + 0.10);
    
            I = d.RT>(ts(i)+0.15) & abs(pulse_size)~=0 & ~isnan(nanmean(datamp(tind,:))');
    
            h = zscore_bygroup(nanmean(datamp(tind,I))', [coh(I)]);
            depvar = h;


            B(i,cont) = nanmean(depvar(pulse_size(I)<0)) - nanmean(depvar(pulse_size(I)>0));
 

        end
    end
end

save('./output_data/beta_pulses','ts','B');





