%analyze 1st batch (9 subjects)
% need Alex's lib
strDirBehav = '/home/liwei/Documents/studies/att_conj/behav';
nSubj = 9;
nSessions = 10;

cPathSessions = cell(nSubj, nSessions);
for kSubj = 1:nSubj
    for kSession = 1:nSessions
        re = ['Training-',num2str(kSubj),'.+s',num2str(kSession),'$'];
        cPathSessions(kSubj, kSession) = FindFiles(strDirBehav, re);
    end
end

cRes = cellfun(@(p) ana_session(p), cPathSessions, 'uni', 0);

sRes = cell2mat(cRes);

t = 1:nSessions;
hl = NaN(nSubj, nSessions);
ll = NaN(nSubj, nSessions);
hlrsvp = NaN(nSubj, nSessions);
llrsvp = NaN(nSubj, nSessions);

for subj = 1:nSubj
    hl(subj,:) = cat(2, sRes(subj,:).hl);
    ll(subj,:) = cat(2, sRes(subj,:).ll);
    hlrsvp(subj,:) = cat(2, sRes(subj,:).hlrsvp);
    llrsvp(subj,:) = cat(2, sRes(subj,:).llrsvp);
end


x = bsxfun(@minus,x,mean(x));

%% plot raw task performance
figure()
plot(t,[hl;ll]);
hold on

phl = plot(t, hl, '^');
pll = plot(t, ll, 'o');
legend([phl(1), pll(1)], 'High', 'Low')


title('Learning')
ylabel('Accuracy')
xlabel('Session')
hold off


%% plot RSVP
figure()
plot(t,[hlrsvp;llrsvp]);
hold on

phlrsvp = plot(t, hlrsvp, '^');
pllrsvp = plot(t, llrsvp, 'o');
legend([phlrsvp(1), pllrsvp(1)], 'HighRSVP', 'LowRSVP')


title('RSVP')
ylabel('Accuracy')
xlabel('Session')
hold off

%% correlation analysis
hlcor = corrcoef([hl',(1:10)']);
hlcor = hlcor(1:(end-1),end);
llcor = corrcoef([ll',(1:10)']);
llcor = llcor(1:(end-1),end);
[H,P,CI,STATS] = ttest(hlcor-llcor);

hlrsvpcor = corrcoef([hlrsvp',(1:10)']);
hlrsvpcor = hlrsvpcor(1:(end-1),end);
llrsvpcor = corrcoef([llrsvp',(1:10)']);
llrsvpcor = llrsvpcor(1:(end-1),end);
[H,P,CI,STATS] = ttest(hlrsvpcor-llrsvpcor);
% 
% H =
% 
%      0
% 
% 
% P =
% 
%     0.5848
% 
% 
% CI =
% 
%    -0.2402
%     0.1451
% 
% 
% STATS = 
% 
%     tstat: -0.5693
%        df: 8
%        sd: 0.2506

%% LME

Y = [hl; ll];
X = cell(18,2);
X(1:8,1) = {'High'};
X(9:end,1) = {'Low'};

%% anova to validate control

% diffrsvp = (llrsvp-hlrsvp)';
% anova1(diffrsvp);
% 
% anova2([hlrsvp;llrsvp],9);

% rsvp = [hlrsvp;llrsvp]';
% loads = [repmat({'hl'},size(hlrsvp,1),1); repmat({'ll'},size(llrsvp,1),1)]';
% subjs = repmat(num2cell('1':'9'),1,2);
% group = [loads;subjs];
% anova1(rsvp, loads)


rsvp = reshape([hlrsvp;llrsvp]',[],1);
loads = reshape([repmat({'hl'},size(hlrsvp)); repmat({'ll'},size(llrsvp))]',[],1);
subjs = reshape(repmat(num2cell('1':'9'),10,2),[],1);
sessions = reshape(repmat(num2cell(('0':'9')'),1,18),[],1);
[~,~,stats] = anovan(rsvp,{loads,subjs,sessions});
c = multcompare(stats,'Dimension',1);
