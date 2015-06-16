function [ res ] = ana_session_control( pathdata )
%ANA_SESSION_CONTROL Summary of this function goes here
%   Detailed explanation goes here
data = readtable(pathdata, 'delimiter', 'tab');

datahlmask = ismember(data.load, 'h') & data.block < 7;
datallmask = ismember(data.load, 'l') & data.block < 7;

res.hlrsvp = mean(data.cor(datahlmask));
res.llrsvp = mean(data.cor(datallmask));
res.hduration = mean(data.duration(datahlmask));
res.lduration = mean(data.duration(datallmask));

end

