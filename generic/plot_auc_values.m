function [ind, x, m, stde, labels] = plot_auc_values(auc, auc_stde, L, str, varargin)

ind = [];
for i=1:length(varargin)
    if isequal(varargin{i},'sort')
        ind = varargin{i+1};
    end
end

for i=1:length(L)
    I(i) = find(ismember(str,L{i}));
end

    
auc = auc(I);

if isempty(ind)
    % [~,ind] = sort(auc,'ascend');
    ind = 1:length(L);
end

m = auc(ind);
stde = auc_stde(I);

x = 1:length(ind);
labels = L(ind);


hold all

terrorbar(x, m, stde,'marker','o','markerfacecolor',0.7*[1,1,1],'color','k','LineStyle','none');
grid on
set(gca,'xticklabel',labels,'xtick',1:length(ind));

xlim([0.8, length(ind)+0.2]);



end