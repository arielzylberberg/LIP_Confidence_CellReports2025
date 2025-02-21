addpath('../generic');

%%

load ../03_Kmeans_calculations/output_data/out_3clusters.mat
load('./output_data/beta_pulses');


%% one approach: bootstrapping

rng(4414,'twister');
nboot = 5000;
for i=1:3
    J = find(clusters==i);

    for k=1:nboot
        n = length(J);
        I = randsample(J,n,1);

        y = nanmean(B(:,I),2);

        % normalize the weights
        tind = ts>=0 & ts<=0.2;
        y = y - nanmean(y(tind));


        [c2,gof,vals(i),fallo,fit_parms_aux(i,:)] = fit_PRR_function(ts,y);
        fit_alpha_boot(i, k) = fit_parms_aux(i,1);
    end

end


%%
save('./output_data/out_boot_shuffle','fit_alpha_boot');

