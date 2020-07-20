clear

global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;


fundir = [pwd];
datadir = [fundir '/../data/PEST/'];
PHITdir = [fundir '/../data/PHIT/'];

%% load PEST data
d = dir([datadir '*.dat']);
for sn = 1:length(d)
    [head, sub(sn)] = readFileSubject_v2([datadir d(sn).name], 'young');
end
L = length(sub);

%% load PHIT data
phit = loadPHIT_v1([PHITdir 'US1_Final_Activity_21days.xlsx']);

%% get unique PHIT emails
EM = vertcat(phit.emails);
phitEmails = unique(EM);

%% demographics for PEST
disp(' ')
disp('demographics for PEST =============================================')
% number of women
disp(['n women = ' num2str(sum(strcmp({sub.sex}', 'F')))])
disp(['n men = ' num2str(sum(strcmp({sub.sex}', 'M')))])

% mean age
disp(['mean age = ' num2str(mean([sub.age]))]);
disp(['max age = ' num2str(max([sub.age]))]);
disp(['min age = ' num2str(min([sub.age]))]);
disp(['std age = ' num2str(std([sub.age]))]);
disp('===================================================================')
disp(' ')

%% ========================================================================
%% BASIC BEHAVIOR ON PEST %% BASIC BEHAVIOR ON PEST %%
%% BASIC BEHAVIOR ON PEST %% BASIC BEHAVIOR ON PEST %%
%% BASIC BEHAVIOR ON PEST %% BASIC BEHAVIOR ON PEST %%
%% ========================================================================

%% BASIC ACCURACY
for sn = 1:length(sub)
    
    safeId = strcmp(sub(sn).realId, 'Safe');
    scamId = strcmp(sub(sn).realId, 'Scam');
    
    realId = strcmp(sub(sn).type, 'Type');
    fakeId = ~realId;
    
    % false positives (safe classified as scam)
    safeScore(sn) = nanmean(sub(sn).userId(safeId));
    FP(sn) = nanmean(sub(sn).userId(safeId) > 2.5);
    
    % false negatives (scam classified as safe)
    scamScore(sn) = nanmean(sub(sn).userId(scamId));
    FN(sn) = nanmean(sub(sn).userId(scamId) < 2.5);
    
    % overall accuracy
    ACC(sn) = nanmean([
        (sub(sn).userId(scamId) > 2.5);
        (sub(sn).userId(safeId) < 2.5)]);
    
    % accuracy by type for the four types
    ACC_safeReal(sn) = nanmean(sub(sn).userId(safeId&realId) < 2.5);
    ACC_safeFake(sn) = nanmean(sub(sn).userId(safeId&fakeId) < 2.5);
    ACC_scamReal(sn) = nanmean(sub(sn).userId(scamId&realId) > 2.5);
    ACC_scamFake(sn) = nanmean(sub(sn).userId(scamId&fakeId) > 2.5);
    
end


AA = [ACC_safeReal' ACC_safeFake' ACC_scamReal' ACC_scamFake'];

%% Figure 2 - Suspicion scores in PEST
clear score
for sn = 1:length(sub)
    iReal = strcmp(sub(sn).type, 'Type');
    iFake = ~iReal;
    iScam = strcmp(sub(sn).realId, 'Scam');
    iSafe = strcmp(sub(sn).realId, 'Safe');
    
    scoreFakeSafe(sn) = nanmean(sub(sn).userId(iFake & iSafe));
    scoreFakeScam(sn) = nanmean(sub(sn).userId(iFake & iScam));
    scoreRealSafe(sn) = nanmean(sub(sn).userId(iReal & iSafe));
    scoreRealScam(sn) = nanmean(sub(sn).userId(iReal & iScam));
    
    score(sn,1) = nanmean(sub(sn).userId(iFake & iSafe));
    score(sn,2) = nanmean(sub(sn).userId(iFake & iScam));
    score(sn,3) = nanmean(sub(sn).userId(iReal & iSafe));
    score(sn,4) = nanmean(sub(sn).userId(iReal & iScam));
end

lab = {'simulated-safe' 'simulated-phish' 'real-safe' 'real-phish'};

order = [3 1 4 2 ];

figure(1); clf;
set(gcf, 'Position', [440   378   580   550])
hg = [0.52 0.15];
wg = [0.15  0.03];
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);
ax(2:3) = easy_gridOfEqualFigures([0.12 0.57], [0.12 0.15 0.03]);

set(ax, 'fontsize', 18, 'tickdir', 'out')
clear l
axes(ax(1)); hold on;
for i = 1:4
    l(i) = plot(i+(rand(size(score,1),1)-0.5)/4,score(:,order(i)),'o', 'markersize', 3);
end

[b] = boxplot(score(:,order), 'notch', 'on', 'outliers', 1, 'widths', 0.4);
set(b, 'linewidth', 2)
set(b(7,:), 'visible','off')
set(b(:,1), 'color', AZblue)
set(b(:,2), 'color', AZblue)
set(b(:,3), 'color', AZred)
set(b(:,4), 'color', AZred)

set(l(1:2), 'color', AZblue*0.5+0.5)
set(l(3:4), 'color', AZred*0.5+0.5)

ylim([1 5])
set(gca, 'view', [90 90],'box', 'off', ...
    'xtick', [1:4], 'xticklabel', lab(order), 'ytick', 1:4, ...
    'yaxislocation', 'right');


mm = max(score(:,order));

plot( [1 1 2 2], [mm(1)+0.1 4.2 4.2 mm(2)+0.1], 'k-')
plot( [3 3 4 4], [mm(3)+0.1 4.2 4.2 mm(4)+0.1], 'k-')
plot( [1.5 1.5 3.5 3.5], [4.5 4.9 4.9 4.7], 'k-')
text(1.7, 4.25, '*', 'fontsize', 30)
text(3.45, 4.25, 'n.s.', 'fontsize', 18)
text(2.7, 4.95, '*', 'fontsize', 30)


dum = get(ax(1), 'position');
wg(1) = dum(1);
wb = dum(3);
annotation('textbox', [wg(1) hg(1)+hb(1)+0.05 wb(1)/4*3 hg(2)], 'string', 'mean suspicion score', 'horizontalalignment', 'center', 'fontsize', 22, 'linestyle', 'none', 'verticalalignment', 'middle')

dw = 0.1;
dh = 0.14;
annotation('textbox', [wg(1)-dw hg(1)+hb(1) 2*dw dh], 'string', sprintf('definitely\nsafe'), 'horizontalalignment', 'center', 'fontsize', 18, 'linestyle', 'none', 'verticalalignment', 'middle')
annotation('textbox', [wg(1)-dw+wb(1)/4*3 hg(1)+hb(1) 2*dw dh], 'string', sprintf('definitely\nsuspicious'), 'horizontalalignment', 'center', 'fontsize', 18, 'linestyle', 'none', 'verticalalignment', 'middle')

dw = 0.1
a = annotation('doublearrow', [wg(1)+dw wg(1)+wb(1)/4*3-dw], [hg(1)+hb(1)+dh/2]*[1 1], 'linewidth', 5, 'head1length', 20, 'head1width', 20, 'head2length', 20, 'head2width', 20);


% fake-safe vs real-safe
axes(ax(2)); hold on;
i = 1; j = 3;
plot(score(:,i), score(:,j), 'o', 'markersize', 7, 'color', AZblue);%*0.5+0.5)
plot([1 4], [1 4], 'k--')
xlim([1 4]); ylim([1 4])
xlabel(lab{i});
ylabel(lab{j})
ll = lsline;
set(ll, 'linewidth', 6, 'color', AZblue)
[r,p] = nancorr(score(:,i), score(:,j), 'spearman');
text(1.1, 3.8, sprintf('r(%d) = %.2f', size(score,1)-2, r),...
    'fontsize', 16)

% real-safe vs fake-safe
axes(ax(3)); hold on;
i = 2; j = 4;
plot(score(:,i), score(:,j), 'o', 'markersize', 7, 'color', AZred);
plot([1 4], [1 4], 'k--')
xlim([1 4]); ylim([1 4])
xlabel(lab{i});
ylabel(lab{j})
ll = lsline;
set(ll, 'linewidth', 6, 'color', AZred)
[r,p] = nancorr(score(:,i), score(:,j), 'spearman');
text(1.1, 3.8, sprintf('r(%d) = %.2f', size(score,1)-2, r),...
    'fontsize', 16)

addABCs(ax(1), [-0.21 0.03], 28)
addABCs(ax(2:3), [-0.08 0.04], 28, 'BC')
saveFigurePdf(gcf, '~/Desktop/Figure_2')
saveFigurePng(gcf, '~/Desktop/Figure_2')

%% ANOVA

S = repmat([1:size(score,1)]', [1 size(score,2)]);
% legitimacy - phish/safe
L = [ones(size(score,1),1) ones(size(score,1),1)*2 ones(size(score,1),1) ones(size(score,1),1)*2];
% authorship - real/constructed
A = [ones(size(score,1),1) ones(size(score,1),1) ones(size(score,1),1)*2 ones(size(score,1),1)*2];

[p,t,stats] = anovan(score(:), [L(:) A(:) S(:)], 'random', 3, ...
    'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0; ], ...
    'varnames', {'legit' 'source' 'subject'});

%%
iError = find(strcmp({t{:,1}}, 'Error'));
iLegit = find(strcmp({t{:,1}}, 'legit'));
iSouce = find(strcmp({t{:,1}}, 'source'));
iInt = find(strcmp({t{:,1}}, 'legit*source'));

partialEta_legit = t{iLegit,2} / (t{iLegit,2} + t{iError,2})
partialEta_author = t{iSouce,2} / (t{iSouce,2} + t{iError,2})
partialEta_Int = t{iInt,2} / (t{iInt,2} + t{iError,2})

%% ANOVA WITH GENDER

S = repmat([1:size(score,1)]', [1 size(score,2)]);
% legitimacy - phish/safe
L = [ones(size(score,1),1) ones(size(score,1),1)*2 ones(size(score,1),1) ones(size(score,1),1)*2];
% authorship - real/constructed
A = [ones(size(score,1),1) ones(size(score,1),1) ones(size(score,1),1)*2 ones(size(score,1),1)*2];
% gender
G = 1+ strcmp({sub.sex}', 'F');
G = [G G G G];

anovan(score(:), [L(:) A(:) G(:) S(:)], 'random', 4, ...
    'nested', [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 1 0], ...
    'model', [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 0 0; 1 0 1 0; 0 1 1 0; 1 1 1 0], ...
    'varnames', {'legit' 'source' 'gender' 'subject'});

%% posthoc t-tests
% post-hoc t-tests
disp('real-safe vs simulated-safe');%: t(96) = XXX, p = XXX,
[~,p,~,stats] = ttest(score(:,1), score(:,3))

disp('real-phish vs simulated-phish');%: t = XXX, p = XXX
[~,p,~,stats] = ttest(score(:,2), score(:,4))

%% accuracy difference between men and women
[~,p, ~, stats] = ttest2(ACC(G(:,1)==1),ACC(G(:,1)==2))
nanmean(ACC(G(:,1)==1))
nanmean(ACC(G(:,1)==2))

%% Supplementary Figure 2 - Mean suspicion score as a function of gender
%% (male vs. female) and email type (real-safe, simulated-safe, real-phish,
%% simulated-phish)
G = strcmp({sub.sex}', 'F');

F = score(G,:);
M = score(~G,:);


clear l
MM = [M; nan(size(F,1)-size(M,1), 4)]
YY = [MM(:,order) F(:,order)];
ZZ = YY(:, [1 5 2 6 3 7 4 8]);
figure(1); clf;
set(gcf, 'Position', [440   378   550   400])
hg = [0.03 0.25];
wg = [0.2 0.09];
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);

hold on;

for i = 1:size(ZZ,2)
    l(i) = plot(i + (rand(size(ZZ,1),1)-0.5)/4, ZZ(:,i),'o');
end
b = boxplot(ZZ, 'notch', 'on');

set(l(1:4), 'color', AZblue*0.5+0.5)
set(l(5:8), 'color', AZred*0.5+0.5)

set(gca, 'view', [90 90], 'xtick', [1:8], ...
    'xticklabel', {
    'male-real-safe'
    'female-real-safe'
    'male-simulated-safe'
    'female-simulated-safe'
    'male-real-phish'
    'female-real-phish'
    'male-simulated-phish'
    'female-simulated-phish'
    }, 'tickdir', 'out', 'box', 'off', 'fontsize', 18, ...
    'yaxislocation', 'right')
dum = get(ax(1), 'position');
wg(1) = dum(1);
wb = dum(3);
annotation('textbox', [wg(1) hg(1)+hb(1)+0.1 wb(1) hg(2)], 'string', 'mean suspicion score', 'horizontalalignment', 'center', 'fontsize', 22, 'linestyle', 'none', 'verticalalignment', 'middle')

dw = 0.1;
dh = 0.25;
annotation('textbox', [wg(1)-dw hg(1)+hb(1) 2*dw dh], 'string', sprintf('definitely\nsafe'), 'horizontalalignment', 'center', 'fontsize', 18, 'linestyle', 'none', 'verticalalignment', 'middle')
annotation('textbox', [wg(1)-dw+wb(1) hg(1)+hb(1) 2*dw dh], 'string', sprintf('definitely\nsuspicious'), 'horizontalalignment', 'center', 'fontsize', 18, 'linestyle', 'none', 'verticalalignment', 'middle')

dw = 0.1
a = annotation('doublearrow', [wg(1)+dw wg(1)+wb(1)-dw], [hg(1)+hb(1)+dh/2]*[1 1], 'linewidth', 5, 'head1length', 20, 'head1width', 20, 'head2length', 20, 'head2width', 20);

ylim([1 4])
plot([0 9],[1 1]*2.5,  'k--')
set(b, 'linewidth', 2)
set(b(7,:), 'visible','off')
set(b(:,1:4), 'color', AZblue)
set(b(:,5:8), 'color', AZred)

saveFigurePdf(gcf, '~/Desktop/SupplementaryFigure_2')
saveFigurePng(gcf, '~/Desktop/SupplementaryFigure_2')

%% ttests
% men showing increased suspicion scores for safe emails (STATS) and
idx = [1 3];
[~,p] = ttest2(nanmean(F(:,idx),2), nanmean(M(:,idx),2))



%% breakout emails by type - constructed/genuine and safe/phish
clear score
for sn = 1:length(sub)
    iReal = strcmp(sub(sn).type, 'Type');
    iFake = ~iReal;
    iScam = strcmp(sub(sn).realId, 'Scam');
    iSafe = strcmp(sub(sn).realId, 'Safe');
    
    scoreFakeSafe(sn) = nanmean(sub(sn).userId(iFake & iSafe));
    scoreFakeScam(sn) = nanmean(sub(sn).userId(iFake & iScam));
    scoreRealSafe(sn) = nanmean(sub(sn).userId(iReal & iSafe));
    scoreRealScam(sn) = nanmean(sub(sn).userId(iReal & iScam));
    
    score(sn,1) = nanmean(sub(sn).userId(iFake & iSafe));
    score(sn,2) = nanmean(sub(sn).userId(iFake & iScam));
    score(sn,3) = nanmean(sub(sn).userId(iReal & iSafe));
    score(sn,4) = nanmean(sub(sn).userId(iReal & iScam));
end

lab = {'fake-safe' 'fake-scam' 'real-safe' 'real-scam'};


%% Figure 3 - Average suspicion scores for each email type
EM = vertcat(sub.emailCode);
TY = strcmp(vertcat(sub.type), 'Type'); % 1 for real
ID = strcmp(vertcat(sub.realId), 'Scam'); % 1 for scam
[sitEmails, i1, i2] = unique(EM);
ty = TY(i1); id = ID(i1);
[sitScore, sitN, sitRT] = scoreSIT_v1(sub, sitEmails);

S_safeReal = (sitScore( (id==0) & (ty==1) ));
S_safeFake = (sitScore( (id==0) & (ty==0) ));
S_scamReal = (sitScore( (id==1) & (ty==1) ));
S_scamFake = (sitScore( (id==1) & (ty==0) ));

S = [sort(S_safeReal); sort(S_safeFake); sort(S_scamReal); sort(S_scamFake)];
I = [ones(size(S_safeReal))
    ones(size(S_safeFake))*2
    ones(size(S_scamReal))*3
    ones(size(S_scamFake))*4]
X = 1:length(S);

clear l xt
figure(1); clf;
set(gcf, 'Position', [440   306   500   600])
hg = [0.03 0.18];
wg = [0.17 0.1];
[ax, hb, wb] = easy_gridOfEqualFigures(hg, wg);
hold on;
for i = 1:4
    l(i) = plot(X(I==i), S(I==i),'.');
    xt(i) = max(X(I==i))+1;
end
xt = [1 xt];
xtt = xt;
xt(end) = xt(end)-1;
set(l(1:2), 'color', AZblue)
set(l(3:4), 'color', AZred)
set(l, 'markersize', 30)
for i = 2:5
    plot(xt(i)*[1 1], [1 4], 'k--')
end

set(gca, 'view', [90 90], 'tickdir', 'out', ...
    'xtick', xt, 'fontsize', 18, 'ytick', [1:4], ...
    'yaxislocation', 'right')

annotation('textbox', [wg(1) hg(1)+hb(1)+0.05 wb(1) hg(2)], 'string', 'mean suspicion score', 'horizontalalignment', 'center', 'fontsize', 22, 'linestyle', 'none', 'verticalalignment', 'middle')
xlabel('email number', 'fontsize', 22)
xlim([0 xt(end)+1])
t = text(mean(xt(1))+15, 4, 'real-safe');
t(2) = text(mean(xt(2))+15, 4, 'simulated-safe');
t(3) = text(mean(xt(3))+15, 4, 'real-phish');
t(4) = text(mean(xtt(4))+15, 4, 'simulated-phish');
set(t, 'fontsize', 18, 'horizontalalignment', 'right')

dw = 0.1;
dh = 0.18
annotation('textbox', [wg(1)-dw hg(1)+hb(1) 2*dw dh], 'string', sprintf('definitely\nsafe'), 'horizontalalignment', 'center', 'fontsize', 18, 'linestyle', 'none', 'verticalalignment', 'middle')
annotation('textbox', [wg(1)-dw+wb(1) hg(1)+hb(1) 2*dw dh], 'string', sprintf('definitely\nsuspicious'), 'horizontalalignment', 'center', 'fontsize', 18, 'linestyle', 'none', 'verticalalignment', 'middle')

dw = 0.1
a = annotation('doublearrow', [wg(1)+dw wg(1)+wb(1)-dw], [hg(1)+hb(1)+dh/2]*[1 1], 'linewidth', 5, 'head1length', 20, 'head1width', 20, 'head2length', 20, 'head2width', 20);
saveFigurePng(gcf, '~/Desktop/Figure_3')
saveFigurePdf(gcf, '~/Desktop/Figure_3')

%% how many emails of each type?

N_safeReal = sum((id==0) & (ty==1) )
N_safeFake = sum((id==0) & (ty==0) )
N_scamReal = sum((id==1) & (ty==1) )
N_scamFake = sum((id==1) & (ty==0) )

%% stick score onto each email
EM = vertcat(sub.emailCode);
TY = strcmp(vertcat(sub.type), 'Type'); % 1 for real
ID = strcmp(vertcat(sub.realId), 'Scam'); % 1 for scam
[sitEmails, i1, i2] = unique(EM);
ty = TY(i1); id = ID(i1);
[sitScore, sitN, sitRT] = scoreSIT_v1(sub, sitEmails);

for sn = 1:length(sub)
    for t = 1:length(sub(sn).emailCode)
        ind = find(strcmp(sitEmails, sub(sn).emailCode{t}));
        sub(sn).score(t,1) = 2*(sitScore(ind)-min(sitScore))/(max(sitScore)-min(sitScore))-1;
        % score is scaled to be between -1 and +1
    end
end
%% how many people got each email?
figure(1); clf;
set(gcf, 'position', [342   514   650   300])
ax = easy_gridOfEqualFigures([0.2 0.1], [0.15 0.15 0.03]);
l = plot_peoplePerEmail_v1(ax, sub)
[l1, l2, l3] = plot_suspicionScoreDistribution_v1(ax(2), sub)

axes(ax(1));
leg = legend({'real-safe' 'simulated-safe' 'real-phish' 'simulated-phish'});
set(leg, 'position', [0.2726    0.6742    0.18    0.2717])
ylabel({'number of' 'participants rating'})
set(ax(1), 'xtick', [1 100:100:200 348])
addABCs(ax, [-0.11 0.09], 32)

saveFigurePdf(gcf, '~/Desktop/SupplementaryFigure_XXX')
saveFigurePng(gcf, '~/Desktop/SupplementaryFigure_XXX')


%% Supplementary Figure 1 ? Number of emails and distribution of suspicion
%% scores for presented emails seen by each participant in PEST.
clear ss
for sn = 1:length(sub)
    [h,x] = hist((sub(sn).score+1)*1.5+1, (1+[-1:0.2:1])*1.5+1);
    ss(sn) = mean(sub(sn).score);
    sub(sn).hScores = h/sum(h);
end


YY = vertcat(sub.hScores);
XX = repmat(x, size(YY,1), 1);
RR = (rand(size(XX))-0.5)/10;

figure(1); clf;
hold on;
plot(RR'+XX',YY', '-','color', [1 1 1]*0.75)
ff = 0.5;
plot(RR'+XX',YY', 'o','color', AZblue*ff+(1-ff), 'markersize', 5)
plot(x,nanmean(vertcat(sub.hScores)), '-','color', AZblue, 'linewidth', 6)
xlim([0.75 4.25])

set(gca, 'xtick', [1:4], 'tickdir', 'out', 'fontsize', 18)
xlabel('suspicion score')
ylabel('frequency')

saveFigurePdf(gcf, '~/Desktop/SupplementaryFigure_1')
saveFigurePng(gcf, '~/Desktop/SupplementaryFigure_1')

%% ========================================================================
%% REGRESSION MODEL %% REGRESSION MODEL %% REGRESSION MODEL %%
%% REGRESSION MODEL %% REGRESSION MODEL %% REGRESSION MODEL %%
%% REGRESSION MODEL %% REGRESSION MODEL %% REGRESSION MODEL %%
%% ========================================================================

%% regression model with stimulus value given by average of everyone's score

%
%
%% Figure 5 - Regression model of PEST behavior
clear B
for sn = 1:length(sub)
    Y = (sub(sn).userId-2.5)/1.5; % scale to be between -1 and +1
    f = sub(sn).score;
    fPast = [nan; f(1:end-1)];
    cPast = [nan; Y(1:end-1)];
    
    
    [B(sn,:)] = glmfit([f fPast cPast], Y);
end

mm = max(B);
mm(3) = min(B(:,3));
figure(1); clf;
set(gcf, 'Position', [440   569   500   220]);
ax = easy_gridOfEqualFigures([0.05 0.1], [0.15 0.03])
hold on;

clear l
for i = 1:size(B,2)
    l(i) = plot(i+(rand(size(B,1),1)-0.5)/4, B(:,i), 'o', 'markersize', 3)
end
b = boxplot(B, 'notch', 'on', 'widths', 0.4);
set(b(7,:), 'visible', 'off')
set(b(:,1:4), 'color', 'k')

set(l(1), 'color', AZred*0.5+0.5)
set(l(2), 'color', AZblue*0.5+0.5)
set(l(3), 'color', AZsand*0.5+0.5)
set(l(4), 'color', AZcactus*0.5+0.5)
set(b, 'linewidth', 2)

text(1.18, mm(1)+0.05, '*', 'fontsize', 36)
text(2.18, mm(2)+0.05, '*', 'fontsize', 36)
text(3.18, mm(3)-0.15, '*', 'fontsize', 36)
text(4.18, mm(4)+0.05, '*', 'fontsize', 36)


plot([0 size(B,2)+1], [0 0], 'k--')
set(gca, 'view', [90 90], 'xtick', 1:4, ...
    'xticklabel', {'phish bias' 'current stimulus' 'past stimulus' 'past rating'}, ...
    'fontsize', 18, 'yaxislocation', 'right', 'tickdir', 'out', ...
    'box', 'off', 'ytick', [-2:0.5:2])
ylim([-0.6 1.6])
ylabel('regression weight')
[~,p, ~, stats] = ttest(B(:,:));

for i = 1:4
    sprintf('t(%d) = %.2f, p = %.2e', stats.df(i), stats.tstat(i), p(i))
end
saveFigurePdf(gcf, '~/Desktop/Figure_5')
saveFigurePng(gcf, '~/Desktop/Figure_5')


%% Supplementary Figure 3 - Regression analysis including terms up to 4
%% trials into the past.
clear B
nPast = 4;
for sn = 1:length(sub)
    clear fPast cPast
    Y = (sub(sn).userId-2.5)/1.5; % scale to be between -1 and +1
    f = sub(sn).score;
    dum = f;
    for n = 1:nPast
        fPast(:,n) = [nan; dum(1:end-1)];
        dum = [nan; dum(1:end-1)];
    end
    dum = Y;
    for n = 1:nPast
        cPast(:,n) = [nan; dum(1:end-1)];
        dum = [nan; dum(1:end-1)];
    end
    
    [B(sn,:)] = glmfit([f fPast cPast], Y);
end

mm = max(B);

figure(1); clf;
set(gcf, 'Position', [440   569   400   400]);

ax =easy_gridOfEqualFigures([0.02 0.1], [0.12 0.05]);

clear l
axes(ax(1)); hold on;
M = nanmean(B);
S = nanstd(B)/sqrt(length(sub));

for i = 1:size(B,2)
    l(i) = plot(i+(rand(size(B,1),1)-0.5)/4, B(:,i), 'o', 'markersize', 3)
end
b = boxplot(B, 'notch', 'on', 'widths', 0.5);
set(b(7,:), 'visible', 'off')
set(b(:,:), 'color', 'k')


set(b(7,:), 'visible', 'off')
set(b(:,:), 'color', [1 1 1]*0)
set(b(6,:), 'color', 'k')

ylabel('regression weight')

plot([0 12], [0 0], 'k--')

set(ax, 'view', [90 90]);

set(l(1), 'color', AZred*0.5+0.5)
set(l(2), 'color', AZblue*0.5+0.5)

set(l(3:nPast+2), 'color', AZsand*0.5+0.5)
set(l(nPast+3:end), 'color', AZcactus*0.5+0.5)
set(b, 'linewidth', 1)

[h,p] = ttest(B);
for i = 1:length(p)
    if p(i) < 0.05
        if p(i) < 0.001
            text(i, mm(i)+0.05, sprintf('p = %.2e', p(i)))
        else
            text(i, mm(i)+0.05, sprintf('p = %.3f', p(i)))
        end
    end
end

set(ax, 'xtick', 1:12, ...
    'fontsize', 18, 'yaxislocation', 'right', 'tickdir', 'out', ...
    'box', 'off', 'ytick', [-2:0.5:2] ,'xtick', [1:10], ...
    'xticklabel', {'phish bias' 'current stimulus' 'stim t-1' 'stim t-2' 'stim t-3' 'stim t-4' 'rating t-1' 'rating t-2' 'rating t-3' 'rating t-4'})

[~,p, ~, stats] = ttest(B(:,:));

for i = 1:4
    sprintf('t(%d) = %.2f, p = %.2e', stats.df(i), stats.tstat(i), p(i))
end
ylim([-0.5 2])
saveFigurePdf(gcf, '~/Desktop/SupplementaryFigure_3')
saveFigurePng(gcf, '~/Desktop/SupplementaryFigure_3')



%% ========================================================================
%% CORRELATION BETWEEN PEST AND PHIT SCORES %%
%% CORRELATION BETWEEN PEST AND PHIT SCORES %%
%% CORRELATION BETWEEN PEST AND PHIT SCORES %%
%% ========================================================================
%% get unique PHIT emails
EM = vertcat(phit.emails);
phitEmails = unique(EM);

%% get PHIT score for each email
[phitScore, phitN] = scorePHIT_v1(phit, phitEmails);
[phitScore_young, phitN_young] = scorePHIT_v1(phit([phit.age]==1), phitEmails);
[phitScore_old, phitN_old] = scorePHIT_v1(phit([phit.age]==2), phitEmails);

%% go through PHIT emails and get PEST score
[sitScore, sitN] = scoreSITPHIT_v1(sub, phitEmails);

%% FIGURE 4 - correlation between PHIT and SIT
ind = (phitScore < 1)
[r,p] = corr(sitScore(ind), phitScore(ind), 'type', 'spearman')

figure(1); clf;
set(gcf, 'Position', [417   518   600   300])
ax = easy_gridOfEqualFigures([0.27 0.12], [0.18 0.53]);
ax(2) = easy_gridOfEqualFigures([0.15 0.16], [0.57 0.02]);
axes(ax(1)); hold on;
plot( sitScore(ind), phitScore(ind), 'o', 'color', [1 1 1]*0.5, ...
    'markersize', 7, 'linewidth', 1)
xlim([1 4])
l = lsline;
set(l, 'color', AZred, 'linewidth', 6)
set(gca, 'xtick', [1:4], 'ytick', [0:0.1:0.5])
xlabel({'mean suspicion score' 'from PEST'})
ylabel({'real-world efficacy' 'in PHIT' '[fraction clicked]'})
text(1.1,0.27, sprintf('r(82) = %.2f\np = %.3f', r, p), 'fontsize', 18, 'interpreter', 'tex')

i1 = phitScore>0;
i0 = phitScore==0;
[~,p,~,stats] = ttest2(sitScore(i1&ind), sitScore(i0&ind))
%

rng(2)
X1 = sitScore(i1&ind);
X0 = sitScore(i0&ind);
X0 = [X0; nan(size(X1,1)-size(X0,1), 1)];

axes(ax(2)); hold on;

l = plot(1+(rand(size(X0))-0.5)/8, X0,'o', 'markersize', 7);
l(2) = plot(2+(rand(size(X1))-0.5)/8, X1,'o', 'markersize', 7);
set(l(1), 'color', AZblue*0.5+0.5)
set(l(2), 'color', AZred*0.5+0.5)
b = boxplot([ X0 X1], 'notch', 'on');
ylim([1 4])
set(gca, 'xtick', [1:2], 'xticklabel', {  'not clicked' 'clicked'})
set(b(:,1), 'color', AZblue)
set(b(:,2), 'color', AZred);
set(b, 'linewidth', 2)
set(gca, 'ytick', [1:4], 'view', [90 90])
set(ax, 'tickdir', 'out', 'fontsize', 18,'box', 'off')
mm = max([X0 X1]);
plot([1 1 2 2], [mm(1)+0.1 3.9 3.9 mm(2)+0.1], 'k')
text(1.4, 4.05, '*', 'fontsize', 30, 'horizontalalignment', 'center')
ylabel({'mean suspicion score' 'from PEST'})
set(ax(2), 'xdir', 'reverse', 'yaxislocation', 'left')

addABCs(ax(1), [-0.078 0.13], 28)
addABCs(ax(2), [-0.12 0.17], 28, 'B')

saveFigurePng(gcf, '~/Desktop/Figure_4')
saveFigurePdf(gcf, '~/Desktop/Figure_4')


