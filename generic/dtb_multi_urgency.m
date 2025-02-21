function [winner,isWinner,dec_time,BupNew,cev] = dtb_multi_urgency(ev,t,a,d,B0,Breflect,USfunc,varargin)

nRaces = length(ev);
N = size(ev{1},3);

if length(B0)==1 %same bound for every condition
    [Bup,~] = expand_bounds(t,B0,a,d,USfunc);
end

%hack because bound has to be positive for mex file
% mBup = min(Bup) - 1;



%decision times
ntr = size(ev{1},1);
decision_step = nan(ntr,nRaces,N);
cev = cell(nRaces,1);
for i=1:nRaces
    e = ev{i};
    cev{i} = nan(size(e,2),size(e,1));
    for k=1:N
        %bounds
        if length(B0)==nRaces
            [Bup,~] = expand_bounds(t,B0(i),a(i),d(i),USfunc);
        end

        BupNew = Bup(1)*ones(size(Bup));
        
        urgency = Bup(1) - Bup;
        durgency = [0; diff(urgency(:))];


        ee = e(:,:,k)';
        ee(1,:) = ee(1,:) - Breflect; %this might generate trouble with params as they are no longer

        ee = ee + durgency(:);
        
        % indep...
        [cev{i}(:,:,k),decision_step(:,i,k)] = single_timebarrier_cross_rectified(ee, BupNew - Breflect,0); %this is a mex file

        
        %this continues to integrate after the bound, but takes longer...
%         [cev{i}(:,:,k),decision_step(:,i,k)] = single_timebarrier_cross_rectified2(ee, BupNew - Breflect,0); %this is a mex file

        cev{i}(:,:,k) = cev{i}(:,:,k) + Breflect;
    end
end

%winner
dt = t(2)-t(1);
winner = nan(ntr,N);
dec_time = nan(ntr,N);
isWinner = nan(ntr,nRaces,N);
for k=1:N
    [winner(:,k),isWinner(:,:,k),dec_time(:,k)] = winning_race(decision_step(:,:,k) * dt);
end


end


function [Bup,Blo] = expand_bounds(t,B0,a,d,USfunc)

t = t(:);
Bup=B0*(t==t);
s=t>d;

switch USfunc
    case 'Linear'
        Bup(s)=B0-a*(t(s)-d);
    case 'Quadratic'
        Bup(s)=B0-a*(t(s)-d).^2;
    case 'Exponential'
        Bup(s)=B0*exp(-a*(t(s)-d));
    case 'Logistic'
        Bup = B0*1./(1+exp(a*(t-d)));
    case 'Hyperbolic'
        Bup(s) = B0*1./(1+(a*(t(s)-d)));
    case 'Step'
        Bup(s)=B0*1e-3;
    case 'Custom1'
        Bup = B0+d*sin(a*t);
    case {'None','Deadline'}
        % Just use Bup itself
    otherwise
        error('USfunc not recognized')
end

% Bounds stop collapsing when the bounds 
% become less than 0.1% of initial height.
% Bup(Bup<=B0*1e-3, 1) = B0*1e-3;




Blo=-Bup;

end


