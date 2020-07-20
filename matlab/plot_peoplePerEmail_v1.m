function l = plot_peoplePerEmail_v1(ax, sub)

emails = {};
eID = [];
for sn = 1:length(sub)
    iReal = strcmp(sub(sn).type, 'Type');
    iFake = ~iReal;
    iScam = strcmp(sub(sn).realId, 'Scam');
    iSafe = strcmp(sub(sn).realId, 'Safe');
    
    emails = {emails{:} sub(sn).emailCode{:}};
    eID = [eID 2*iScam' + iFake' + 1];
end
[uEmails,i,j] = unique(emails);
uID = eID(i);


for sn = 1:length(sub)
    [~,i,j] = intersect(sub(sn).emailCode, uEmails);
    sub(sn).emailNums = j;
    sub(sn).en = zeros(size(uEmails));
    sub(sn).en(j) = 1;
end
XX = vertcat(sub.en);
YY = sum(XX);
clear ZZ
for i = 1:4
    ZZ{i} = YY(uID == i);
end

% figure(1); clf;
axes(ax(1));
hold on;
L = 0;

for i = 1:4
    l(i) = plot([1:length(ZZ{i})]+L, sort(ZZ{i}),'.', 'markersize', 30);
    L = L + length(ZZ{i});
end
xlabel('email number')
ylabel('number of participants rating email')
% legend({'real-safe' 'simulated-safe' 'real-phish' 'simulated-phish'})
ylim([0 100])
xlim([0 L+1])
set(gca, 'tickdir', 'out', 'fontsize', 18, 'xtick', [1 50:50:L L])
