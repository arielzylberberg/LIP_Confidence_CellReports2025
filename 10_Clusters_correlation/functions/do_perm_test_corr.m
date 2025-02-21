function [rho, pval, p_as_extreme] = do_perm_test_corr(X,Y,tind)

% not shuffled
nt = length(tind);
rho = nan(nt);
pval = nan(nt);
for i=1:nt
    for j=1:nt
        x = X(:,tind(i));
        y = Y(:,tind(j));
        K = ~isnan(x) & ~isnan(y);
        [rho(i,j),pval(i,j)] = corr(x(K),y(K));
    end
end


% shuffled
nshuffles = 200;
rho_shuffle = nan(nt,nt,nshuffles);
for m=1:nshuffles
    disp(num2str(m));
    xx = X;
    yy = Y;
    % shuffle rows
    idx = randperm(size(xx,1));
    xx = xx(idx,:);
    for i=1:nt
        for j=1:nt
            x = xx(:,tind(i));
            y = yy(:,tind(j));
            K = ~isnan(x) & ~isnan(y);
            rho_shuffle(i,j,m) = corr(x(K),y(K));
        end
    end
end


% make the stats comparison
% tt = t(tind);
tt = 1:length(tind);
[t1,t2] = meshgrid(tt,tt); %t1: time in ev int dimension; t2: time in momentary ev dimension
Ja = t1(:)<t2(:);
Jb = t1(:)>t2(:);
media = nanmean(rho(Ja))-nanmean(rho(Jb));

% calc on shuffles
media_shuffle = nan(nshuffles,1);
for m=1:nshuffles
    r = rho_shuffle(:,:,m);
    media_shuffle(m) = nanmean(r(Ja))-nanmean(r(Jb));
end
% fit Gaussian to the shuffles
pd = fitdist(media_shuffle,'Normal');
% prob of obtaining a sample as extreme as the one observed
p_as_extreme = 1-pd.cdf(media);


end