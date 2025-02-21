function recolor_psth(pd,color_flag)

if nargin==1
    color_flag = 1;
end
for j=1:length(pd.h_ax)
    h  = get(pd.h_ax(j),'children');

    
    switch color_flag
        case 1
            colores = movshon_colors(6);
            colores = [colores(end:-1:1,:);colores];
            lsty = {'-',':'};
        case 2
            colores = cbrewer('div','BrBG',12);
            lsty = {'-','-'};
        case 3
            colores = cbrewer('div','RdYlBu',12);
            lsty = {'-','-'};
        case 4
            colores = NS_colors(12);
            lsty = {'-','-'};

    end

    for i=1:length(h)
        set(h(i),'color',colores(i,:));
        if i<=6
            set(h(i),'LineStyle',lsty{1});
        else
            set(h(i),'LineStyle',lsty{2});
        end
    end
end