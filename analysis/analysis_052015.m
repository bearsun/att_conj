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


%% plot raw task performance
figure()
phl = plot(t,hl, 'r');
hold on

pll = plot(t,ll, 'g');

legend([phl(1), pll(1)], 'High', 'Low')


title('Learning')
ylabel('Accuracy')
xlabel('Session')
hold off

%% a loop to draw each subj seperatly
figure1 = figure;
for ksubj = 1:nSubj
    subplot(nSubj/3, 3, ksubj);
    plot(t, hl(ksubj, :), 'r');
    hold on
    plot(t, ll(ksubj, :), 'g');
    plot(t, hlrsvp(ksubj, :), 'r--');
    plot(t, llrsvp(ksubj, :), 'g--');
    plot(t, llrsvp(ksubj, :)-hlrsvp(ksubj, :), 'b--');
    plot([0,10], [0 0], 'k');
    
    set(gca, 'Ylim', [-.2, 1]);

    
    title(['Subject ', num2str(ksubj)])
    ylabel('Accuracy')
    xlabel('Session')
    hold off
end

legend(gca, 'HL', 'LL', 'HLrsvp', 'LLrsvp', 'Diffrsvp');

%% center and group average, plot by sessions
aver_hl = mean(bsxfun(@minus, hl, mean(hl,2)), 1);
aver_ll = mean(bsxfun(@minus, ll, mean(ll,2)), 1);
std_hl = std(bsxfun(@minus, hl, mean(hl,2)), 1);
std_ll = std(bsxfun(@minus, ll, mean(ll,2)), 1);

f1 = figure();
p1 = plot(t,[aver_hl; aver_ll]);
hold on
e1 = errorbar(t, aver_hl, std_hl,'LineStyle', 'none', 'Color', p1(1).Color);
e2 = errorbar(t, aver_ll, std_ll,'LineStyle', 'none', 'Color', p1(2).Color);

phl = plot(t, aver_hl, '^', 'Color', p1(1).Color);
pll = plot(t, aver_ll, 'o', 'Color', p1(2).Color);
legend([phl(1), pll(1)], 'High', 'Low')

title('Learning')
ylabel('Accuracy (Centered)')
xlabel('Session')
hold off

%% without subj2
hl_no2 = hl([1,3:end],:);
ll_no2 = ll([1,3:end],:);

aver_hl = mean(bsxfun(@minus, hl_no2, mean(hl_no2,2)), 1);
aver_ll = mean(bsxfun(@minus, ll_no2, mean(ll_no2,2)), 1);
std_hl = std(bsxfun(@minus, hl_no2, mean(hl_no2,2)), 1);
std_ll = std(bsxfun(@minus, ll_no2, mean(ll_no2,2)), 1);

f1 = figure();
p1 = plot(t,[aver_hl; aver_ll]);
hold on
e1 = errorbar(t, aver_hl, std_hl,'LineStyle', 'none', 'Color', p1(1).Color);
e2 = errorbar(t, aver_ll, std_ll,'LineStyle', 'none', 'Color', p1(2).Color);

phl = plot(t, aver_hl, '^', 'Color', p1(1).Color);
pll = plot(t, aver_ll, 'o', 'Color', p1(2).Color);
legend([phl(1), pll(1)], 'High', 'Low')

title('Learning excluding subj2')
ylabel('Accuracy (Centered)')
xlabel('Session')
hold off

%% plot RSVP
figure()
plot(t,hlrsvp, 'Color', 'r');
hold on
plot(t,llrsvp, 'Color', 'g');


phlrsvp = plot(t, hlrsvp, '^', 'Color', 'r');
pllrsvp = plot(t, llrsvp, 'o', 'Color', 'g');
legend([phlrsvp(1), pllrsvp(1)], 'HighRSVP', 'LowRSVP')


title('RSVP')
ylabel('Accuracy')
xlabel('Session')
hold off

%% RSVP diff grouped by subj
drsvp = llrsvp-hlrsvp;
aver_drsvp = mean(drsvp,2);
std_drsvp = std(drsvp,0,2);

figure1 = figure; 
bar1 = bar(aver_drsvp);
hold on
errorbar(1:numel(aver_drsvp), aver_drsvp, std_drsvp, 'LineStyle', 'none');

title('Difference')
ylabel('Accuracy')
xlabel('Subject')
hold off

%% 2nd subject
figure1 = figure;
p2 = plot(1:10, llrsvp(2,:), 'b--');
hold on
p1 = plot(1:10, hlrsvp(2,:),'r--');
p4 = plot(1:10, ll(2,:), 'g');
p3 = plot(1:10, hl(2,:), 'c');

legend([p1(1), p2(1), p3(1), p4(1)], 'HighRSVP', 'LowRSVP', 'High', 'Low')
title('Subject 2')
ylabel('Accuracy')
xlabel('Session')
hold off


%% correlation analysis
hlslope = slope(hl);
llslope = slope(ll);
[H,P,CI,STATS] = ttest(llslopt-hlslope);

[H,P,CI,STATS] = ttest(llslopt);

[H,P,CI,STATS] = ttest(hlslopt);
% hlcor = corrcoef([hl',(1:10)']);
% hlcor = hlcor(1:(end-1),end);
% llcor = corrcoef([ll',(1:10)']);
% llcor = llcor(1:(end-1),end);
% [H,P,CI,STATS] = ttest(hlcor-llcor);
% 
% hlrsvpcor = corrcoef([hlrsvp',(1:10)']);
% hlrsvpcor = hlrsvpcor(1:(end-1),end);
% llrsvpcor = corrcoef([llrsvp',(1:10)']);
% llrsvpcor = llrsvpcor(1:(end-1),end);
% [H,P,CI,STATS] = ttest(hlrsvpcor-llrsvpcor);
% 
% hlcor_no2 = hlcor([1,3:end]);
% llcor_no2 = llcor([1,3:end]);
% [H,P,CI,STATS] = ttest(hlcor_no2-llcor_no2);
% % 
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

%% plot bar
mean_cor = mean([hlslope,llslope])';
std_cor = std([hlslope, llslope])';

figure1 = figure;
axes1 = axes('Parent',figure1,'XTickLabel',{'High Load','Low Load'},...
     'XTick',[1 2]);
box(axes1,'on');
hold(axes1,'on');
 
bar1 = bar(mean_cor);
errorbar(1:2, mean_cor, std_cor, 'LineStyle', 'none');

% mean_cor = [mean([hlcor,llcor]); mean([hlrsvpcor,llrsvpcor])]';
% std_cor = [std([hlcor, llcor]); std([hlrsvpcor, llrsvpcor])]';
% 
% figure1 = figure;
% axes1 = axes('Parent',figure1,'XTickLabel',{'High Load','Low Load'},...
%     'XTick',[1 2]);
% box(axes1,'on');
% hold(axes1,'on');
% 
% bar1 = bar(mean_cor);
% set(bar1(2),'DisplayName','Search');
% set(bar1(1),'DisplayName','RSVP');
% XOffset = cat(2,bar1.XOffset);
% errorbar([1 + XOffset; 2 + XOffset], mean_cor, std_cor, 'LineStyle', 'none');
% 
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.738031914893617 0.719665765642484 0.153739402571097 0.164810768653545],...
%     'FontSize',9);
% 
% bar_cor=bar(mean_cor);
% hold on
% ebar_cor = errorbar(mean_cor, std_cor);
% hold off

%% TL vs RG
tl = NaN(9,10);
rg = NaN(9,10);
tl(1:2:end,:) = hl(1:2:end,:);
tl(2:2:end,:) = ll(2:2:end,:);
rg(1:2:end,:) = ll(1:2:end,:);
rg(2:2:end,:) = hl(2:2:end,:);

figure()
ptl = plot(t, tl, 'r');
hold on
prg = plot(t, rg, 'g');

legend([ptl(1), prg(1)], 'TL', 'RG')


title('Learning')
ylabel('Accuracy')
xlabel('Session')
hold off

aver_tl = mean(bsxfun(@minus, tl, mean(tl,2)), 1);
aver_rg = mean(bsxfun(@minus, rg, mean(rg,2)), 1);
std_tl = std(bsxfun(@minus, tl, mean(tl,2)), 1);
std_rg = std(bsxfun(@minus, rg, mean(rg,2)), 1);

f1 = figure();
p1 = plot(t,[aver_tl; aver_rg]);
hold on
e1 = errorbar(t, aver_tl, std_tl,'LineStyle', 'none', 'Color', p1(1).Color);
e2 = errorbar(t, aver_rg, std_rg,'LineStyle', 'none', 'Color', p1(2).Color);

ptl = plot(t, aver_tl, '^', 'Color', p1(1).Color);
prg = plot(t, aver_rg, 'o', 'Color', p1(2).Color);
legend([ptl(1), prg(1)], 'TL', 'RG')

title('Learning TL vs RG')
ylabel('Accuracy (Centered)')
xlabel('Session')
hold off


tlcor = slope(tl);
rgcor = slope(rg);
[H,P,CI,STATS] = ttest(tlcor,rgcor);

mean_cor = mean([tlcor, rgcor])';
std_cor = std([tlcor, rgcor])';

figure1 = figure;
axes1 = axes('Parent',figure1,'XTickLabel',{'TL','RG'},...
     'XTick',[1 2]);
box(axes1,'on');
hold(axes1,'on');
 
bar1 = bar(mean_cor);
errorbar(1:2, mean_cor, std_cor, 'LineStyle', 'none');
% %% LME
% 
% Y = [hl; ll];
% X = cell(18,2);
% X(1:8,1) = {'High'};
% X(9:end,1) = {'Low'};

%% anova to validate control
% 
% % diffrsvp = (llrsvp-hlrsvp)';
% % anova1(diffrsvp);
% % 
% % anova2([hlrsvp;llrsvp],9);
% 
% % rsvp = [hlrsvp;llrsvp]';
% % loads = [repmat({'hl'},size(hlrsvp,1),1); repmat({'ll'},size(llrsvp,1),1)]';
% % subjs = repmat(num2cell('1':'9'),1,2);
% % group = [loads;subjs];
% % anova1(rsvp, loads)
% 
% 
rsvp = reshape([hlrsvp;llrsvp]',[],1);
loads = reshape([repmat({'hl'},size(hlrsvp)); repmat({'ll'},size(llrsvp))]',[],1);
subjs = reshape(repmat(num2cell('1':'9'),10,2),[],1);
sessions = reshape(repmat(num2cell(('0':'9')'),1,18),[],1);
[~,~,stats] = anovan(rsvp,{loads,subjs,sessions});
c = multcompare(stats,'Dimension',1);
