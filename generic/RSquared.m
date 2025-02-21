function rsquared = RSquared(y,yfit)
% rsquared = RSquared(y,yfit)

TSS = nanmean((y-nanmean(y)).^2);
RSS = nanmean((y-yfit).^2);
rsquared = 1 - RSS/TSS;