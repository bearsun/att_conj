function res = ana_session(pathdata)
% analyze one session of data for att_conj study
%   output res with field hl, ll, hlrsvp, llrsvp, hduration, lduration

data = readtable(pathdata, 'delimiter', 'tab');

datahlmask = ismember(data.load, 'h') & ~isnan(data.rgcor);
datallmask = ismember(data.load, 'l') & ~isnan(data.rgcor);

res.hl = mean(data.rgcor(datahlmask));
res.ll = mean(data.rgcor(datallmask));
res.hlrsvp = mean(data.cor(datahlmask));
res.llrsvp = mean(data.cor(datallmask));
res.hduration = mean(data.duration(datahlmask));
res.lduration = mean(data.duration(datallmask));

end

