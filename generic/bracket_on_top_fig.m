function [hline, hstr] = bracket_on_top_fig(ax_handle, X, multiplier, str)

if nargin==2
    multiplier = 1;
end

pos = get(ax_handle,'Position'); 
ylim = get(ax_handle,'ylim');
xlim = get(ax_handle,'xlim');

delta = 0.1; % distance to top of axis (in % of axis height)

[x_fig1] = ax2fig(ax_handle,X(1),nan);
[x_fig2] = ax2fig(ax_handle,X(2),nan);

y_fig = pos(2) + pos(4);

y_fig = y_fig + pos(4) * delta * multiplier;

hline(1) = annotation('line',[x_fig1, x_fig2],[y_fig,y_fig]);

hstr = annotation('textbox',[x_fig1, y_fig,x_fig2 - x_fig1, 0.05]);

set(hstr,'String',str,'verticalalignment','baseline','horizontalalignment',...
    'center','LineStyle','none','FontAngle','italic');

h = 0.05 * pos(4); % height
hline(2) = annotation('line',[x_fig1, x_fig1],[y_fig,y_fig-h]);
hline(3) = annotation('line',[x_fig2, x_fig2],[y_fig,y_fig-h]);










