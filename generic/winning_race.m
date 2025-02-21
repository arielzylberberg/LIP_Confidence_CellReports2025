function [y,isWinner,val] = winning_race(dt)
% function y = winning_race(dt)
%dt [ntrials x races]: decision time of each race. NaN means
%bound was not reached
%y [ntrials x 1]: id of winner. If two are equal, returns the first
%isWinner: to determine ties.

[val,y]       = min(dt,[],2);
y(isnan(val)) = nan;
isWinner      = bsxfun(@eq,dt,val);

