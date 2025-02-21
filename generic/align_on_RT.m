function [H_rt, t_rt] = align_on_RT(H, RT, t, filt)

if nargin<4
    filt = true(size(RT));
end

H(:,t<0.2,:) = nan; % nan the first 200 ms
E = RT(filt);
times = t;
win = [2,0]; % ignored
nneurons = size(H,1);
ntrials = size(H,3);
for i=1:nneurons
    I = filt;
    h = squeeze(H(i,:,filt));
    [datamp,timeslocked] = eventlockedmatc(h,times,E,win);
    if i==1
        H_rt = nan(nneurons, length(timeslocked), ntrials);
        H_rt(i,:,I) = datamp;
    else
        H_rt(i,:,I) = datamp;
    end
end

t_rt = timeslocked;

end