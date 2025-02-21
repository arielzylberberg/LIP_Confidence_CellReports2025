function [out,p] = logit_fit(scoh,resp,vfilt,doplot,newfig,color)

auto_color = 0;
if nargin<=5 || isempty(color)
    auto_color = 1;
end
    

resp = resp(:);

if nargin<4 || isempty(vfilt)
    vfilt = ones(size(resp));
end

for i=1:size(vfilt,2)
    filt = vfilt(:,i);
    b = glmfit([scoh(filt) ones(sum(filt),1)],resp(filt),'binomial','link','logit','constant','off');
    
    mini = min(scoh(filt));
    maxi = max(scoh(filt));
    xf = linspace(mini,maxi,100)';
    yf = glmval(b,[xf ones(size(xf))],'logit','constant','off');
    
    out(i).xfit = xf;
    out(i).yfit = yf;
    
end


p = [];
if doplot
    for i=1:size(vfilt,2)
        if newfig && i==1
            p = publish_plot(1,1);
        end
        
        if auto_color
            color = get(he,'Color');
        end
        
        [t,x,s] = curva_media(resp,scoh,vfilt(:,i),0);
        he = terrorbar(t,x,s,'marker','.','color',color,'LineStyle','none');
        out(i).hdata = he;
        
        hold all
        hf = plot(out(i).xfit,out(i).yfit,'color',color);
        out(i).hfit = hf;
        if newfig
            p.format('FontSize',15);
        end
    end
    
end

