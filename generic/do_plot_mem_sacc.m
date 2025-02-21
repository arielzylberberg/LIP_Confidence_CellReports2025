function [p, ht] = do_plot_mem_sacc(targ_mem, ineuron, uni_target_pos)

colores = cbrewer('div','RdYlGn',100);
% ineuron = 100; 

p = publish_plot(1,1);

Snorm = targ_mem(:,ineuron) - nanmin(targ_mem(:,ineuron));
Snorm = Snorm/nanmax(Snorm)*100;
if ~all(targ_mem(:,ineuron)==0)
    for i=1:size(uni_target_pos,1)
        %              msize = 1 + Snorm(i);
        color_id = findclose(1:size(colores,1), Snorm(i));
        color = colores(color_id(1),:);
        plot(uni_target_pos(i,1),uni_target_pos(i,2),'o',...
            'markersize',8,'color',color,'markerfacecolor',color);
        hold on
    end
end

str = [num2str(round(nanmax(targ_mem(:,ineuron)) - nanmin(targ_mem(:,ineuron)))),' Hz'];
ht = text(0,0,str,'horizontalalignment','center','verticalalignment','middle');

p.format();