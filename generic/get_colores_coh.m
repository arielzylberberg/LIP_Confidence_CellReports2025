function colores = get_colores_coh(N)

addpath '/Users/arielzy/Dropbox/Data/default_matlab_ariel_files/custom_colormaps/customcolormap'

% J = customcolormap(linspace(0,1,5), {'#3b2026','#7c5288','#d5a6f2','#5e762a','#454541'});

J = customcolormap(linspace(0,1,5),{'#594c91', '#b1a5de' ,'#f1e4fe', '#b6bf85', '#81a34a'});

% colorbar; colormap(J);

n = size(J,1);

idx = 1 + round(linspace(0,1,N)*(n-1));

colores = J(idx,:);

end