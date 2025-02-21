function Hnorm = normalize_neural_resp(H, t, include)

if nargin<3
    include = true(size(H,3),1);
end

tind = t>=-0.2 & t<=0;
h = squeeze(nanmean(H(:,tind,include),2));
base_mean = mean(h,2);
base_std = std(h,[],2);

Hnorm = (H - base_mean)./base_std;

end