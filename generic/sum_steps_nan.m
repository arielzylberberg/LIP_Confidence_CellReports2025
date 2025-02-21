function [Y,ti] = sum_steps_nan(X,step,win,dim)
% function [Y,ti] = sum_steps_nan(X,step,win,dim)
% dada una matrix X, hace nansumas en la dimension dim cada step pasos
% incluyendo "win" columnas
% X es de dimension 1 o 2

if nargin<3
    if isvector(X)
        if iscolumn(X)
            dim = 1;
        else
            dim = 2;
        end
    else
        error('specify dim')
    end
end

J = 1:step:(size(X,dim)-win+1);

for i=1:length(J)
    j = J(i);
    if dim==1
        to = min(j+win-1,size(X,1));
        Y(i,:) = nansum(X(j:to,:),1);
    elseif dim==2
        to = min(j+win-1,size(X,2));
        Y(:,i) = nansum(X(:,j:to),2);
    elseif dim==3
        to = min(j+win-1,size(X,3));
        Y(:,:,i) = nansum(X(:,:,j:to),3);
    end
end

%start of each bin
%dJ = J(2)-J(1);
% ti = J;
ti = J + round(win/2);%centered in the window