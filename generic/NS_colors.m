function colors = NS_colors(N)

if nargin==0
    N = 11;
end
plotPar.colors = {[84 48 5]/255, [140 81 10]/255, [191 129 45]/255, [223 194 125]/255, [246 232 195]/255;...
  [0 60 48]/255, [1 102 94]/255, [53 151 143]/255, [128 205 193]/255, [199 234 229]/255;...
  [142,1,82]/255, [197,27,125]/255, [222,119,174]/255, [241,182,218]/255, [253,224,239]/255};

c1 = cat(1,plotPar.colors{2,:});
c2 = cat(1,plotPar.colors{3,:});

if N==12
    colors = [c1; 0.6*[1,1,1]; 0.8*[1,1,1]; c2(end:-1:1,:)];
elseif N==11
    colors = [c1; 0.7*[1,1,1]; c2(end:-1:1,:)];
elseif N==1
    colors = [0,0,1];
else
    colors = rainbow_colors(N);
end

end