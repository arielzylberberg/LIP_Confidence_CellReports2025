addpath('../functions/');
addpath('../generic/');

%%
idataset = 1;
d = get_data(50,[],idataset);


%%
itr = 200;
s = d.spike_times_ms(itr,:);

dt_ms = 1;
[t,spk] = spikeanalysis.spk_to_hist(s,-200,1000*d.RT(itr),dt_ms);

p = publish_plot(1,1);
set(gcf,'Position',[1306   774   258   418]);
spikeanalysis.plot_raster(spk,t/1000);
hold all
plot([0,0],ylim,'r');
plot([d.RT(itr),d.RT(itr)],ylim,'color',[0,0.5,0]);
xlabel('Time from motion onset [s]');
ylabel('Neuron #');

set(gca,'ytick',[1,50,100,150]);
set(gca,'xtick',[0,.50])
p.offset_axes();
% p.format('FontSize',18);

p.append_to_pdf('fig_raster',1,1);