function [winner,isWinner,dec_time] = dtb_multi(ev,t,a,d,B0,USfunc)

nRaces = length(ev);
N = size(ev{1},3);

if length(B0)==1 %same bound for every condition
    [Bup,~] = expand_bounds(t,B0,a,d,USfunc);
end

%decision times
ntr             = size(ev{1},1);
decision_step   = nan(ntr,nRaces,N);
for i=1:nRaces
    e = ev{i};
    for k=1:N
        %bounds
        if length(B0)==nRaces
            [Bup,~] = expand_bounds(t,B0(i),a(i),d(i),USfunc);
        end
        [~,decision_step(:,i,k)] = single_timebarrier_cross2(e(:,:,k)',Bup,1); %this is a mex file
    end
end

%winner
dt          = t(2)-t(1);
winner      = nan(ntr,N);
dec_time    = nan(ntr,N);
isWinner    = nan(ntr,nRaces,N);
for k = 1:N
    [winner(:,k),isWinner(:,:,k),dec_time(:,k)] = winning_race(decision_step(:,:,k) * dt);
end



end


