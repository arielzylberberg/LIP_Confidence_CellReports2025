function [t,H] = spk_to_hist(spk,tini,tend,dt)
% function H = spk_to_hist(spk,tini,tend,bin_size)


%% 
if nargin<2 
    concatenado = cat(1,spk{:});
    tini = nanmin(concatenado);
    tend = nanmax(concatenado);
    dt = 1;
end

%%
t = tini:dt:tend;
nt = length(t);
ntr = length(spk);
H = zeros(ntr,nt);
for i=1:length(spk)
    if ~isempty(spk{i})
        H(i,:) = histc(spk{i},t);
    end
    
end



end