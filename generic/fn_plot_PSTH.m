function [p, av] = fn_plot_PSTH(t,H,RT,coh,correct)

use_bigSbigT_flag = 1;
if use_bigSbigT_flag

    fn_do_plot_aligned = @(id_neurons_to_include, RT, t, H, coh, I) ...
        (do_plot_aligned_bigSbigT(id_neurons_to_include, RT, t, H, coh, I));
else
    fn_do_plot_aligned = @(id_neurons_to_include, RT, t, H, coh, I) ...
        (do_plot_aligned(id_neurons_to_include, RT, t, H, coh, I));
end

%%



K = correct==1;

[p, av] = fn_do_plot_aligned([], RT(K), t, H(K,:)', coh(K), ones(sum(K),1)==1);

p.current_ax(1);
xlabel({'Time from','motion onset [s]'});

p.current_ax(2);
xlabel({'Time from','response [s]'});

set(p.h_ax(1),'xlim',[-0.2,1]);
set(p.h_ax(2),'xlim',[-0.7,-0.05]);
p.same_ylim_by_row();
p.same_xscale();

% set(p.h_ax(2),'YAxisLocation','right')


%% change colors
colores = movshon_colors(6);
colores = [colores(end:-1:1,:); colores];
lsty = {'-',':'};

for k=1:length(p.h_ax)
    h = get(p.h_ax(k),'children');
    for i=1:length(h)
        set(h(i),'color',colores(i,:));
        if i<=6
            set(h(i),'LineStyle',lsty{1});
        else
            set(h(i),'LineStyle',lsty{2});
        end
    end
end
%%
p.current_ax(1);
hold all
plot([0,0],ylim,'color',0.7*[1,1,1]);



%%
set(p.h_ax(1),'xtick',[0,0.25,0.5,0.75],'xticklabel',[0,0.25,0.5,0.75]);
set(p.h_ax(2),'xtick',[-0.5:0.25:0],'xticklabel',[-0.5:0.25:0],'ycolor','none');



p.format('FontSize',18,'LineWidthPlot',1.0);

end