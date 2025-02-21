function [vals_cat_sym,vals_cat_mean] = symmetric_scale_prctiles(dv,ncats_signed)

if ~iseven(ncats_signed)
    error('ncats has to be even')
end

ncats_unsigned = (ncats_signed)/2;

adv = abs(dv);
[adv_sort,idx_sort] = sort(adv);

ntr = length(dv);
tr_per_cat = ceil(ntr/ncats_unsigned);

bin_idx = repmat([1:ncats_unsigned],tr_per_cat,1);%which trials go to which category, unsigned
bin_idx = bin_idx(1:ntr)';


% bin edges
I1 = [tr_per_cat:tr_per_cat:tr_per_cat*(ncats_unsigned-1)];
I2 = 1 + I1;
edges = [0; [adv_sort(I1) + adv_sort(I2)]/2 ; nanmax(adv)]; 


%unsort
[~,idx_unsort] = sort(idx_sort);
bin_idx = bin_idx(idx_unsort);

% make additional cats for the negative values of dv
bin_idx(dv<0) = -1*bin_idx(dv<0);

% converts to ordered numbers
[~,~,bin_idx] = unique(bin_idx);

bin_x = edges(1:end-1)+diff(edges)/2;
%make the first zero
% bin_x(1) = 0;
%both signs
bin_x = sort([-1*bin_x;bin_x]);

vals_cat_sym = bin_x(bin_idx);

[~,uu] = curva_media(dv,bin_idx,[],0);


vals_cat_mean = index_to_val(bin_idx,uu);

% curva_media(adv,vals,[],1)


