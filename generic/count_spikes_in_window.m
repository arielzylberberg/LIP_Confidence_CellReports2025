function count = count_spikes_in_window(spike_times, idx_neurons, t_from, t_to)

if isempty(idx_neurons)
    idx_neurons = 1:size(spike_times,2);
end

% spkt = D.spike_times_ms(:,idx_neurons);
spkt = spike_times(:,idx_neurons);
[ntrials, nneurons] = size(spkt);

if isscalar(t_from)
    t_from = ones(ntrials,1)*t_from;
end

if isscalar(t_to)
    t_to = ones(ntrials,1)*t_to;
end

count = nan(ntrials,nneurons);
for i=1:ntrials
    for j=1:nneurons
        s = spkt{i,j};
        count(i,j) = sum(s>=t_from(i) & s<=t_to(i));
    end
end

end


