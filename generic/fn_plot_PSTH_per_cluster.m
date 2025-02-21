function [p, av] = fn_plot_PSTH_per_cluster(t,H,RT,coh,correct,clusters, neuron_id)

use_bigSbigT_flag = 1;
if use_bigSbigT_flag

    fn_do_plot_aligned = @(id_neurons_to_include, RT, t, H, coh, I) ...
        (do_plot_aligned_bigSbigT(id_neurons_to_include, RT, t, H, coh, I));
else
    fn_do_plot_aligned = @(id_neurons_to_include, RT, t, H, coh, I) ...
        (do_plot_aligned(id_neurons_to_include, RT, t, H, coh, I));
end

%% 
if numel(clusters)==numel(RT)
    clust = clusters; % already indicated to which cluster each trial belongs to
else
    clust = clusters(neuron_id);
end

n = length(nanunique(clust));
p = publish_plot(n+1,2);
% set(gcf,'Position',[440  180  512  618]);
set(p.h_fig,'Position',[440  180  400  618]);
for i=1:n+1
    I = correct==1;
    if i==1
        K = I;
    else
        K = I & clust==(i-1);
    end

    [p_aux,av(i)] = fn_do_plot_aligned([], RT(K), t, H(K,:)', coh(K), I(K));

    p_aux.current_ax(1);
    if i==1
        title('All neurons');
    else
        nneurons = length(nanunique(neuron_id(K)));
        title(['N = ',num2str(nneurons)]);
    end
    p.copy_from_ax(p_aux.h_ax(1),p.h_ax(2*i-1));
    p.copy_from_ax(p_aux.h_ax(2),p.h_ax(2*i));
    close(p_aux.h_fig);

end



p.current_ax(2*n+1);
xlabel({'Time from','motion onset [s]'});

p.current_ax(2*n+2);
xlabel({'Time from','response [s]'});

set(p.h_ax(1:2:2*(n+1)),'xlim',[-0.2,1]);
set(p.h_ax(2:2:2*(n+1)),'xlim',[-0.7,-0.05]);
p.same_ylim_by_row();
p.unlabel_center_plots();
p.same_xscale();

set(p.h_ax(2:2:end),'YAxisLocation','right')
p.unlabel_center_plots();
p.format('FontSize',12,'LineWidthPlot',1.0);

end