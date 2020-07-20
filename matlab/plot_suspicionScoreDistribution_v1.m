function [l1, l2, l3] = plot_suspicionScoreDistribution_v1(ax, sub)

global AZblue

% clear ss
for sn = 1:length(sub)
    [h,x] = hist((sub(sn).score+1)*1.5+1, (1+[-1:0.2:1])*1.5+1);
    ss(sn) = mean(sub(sn).score);
    sub(sn).hScores = h/sum(h);
end


YY = vertcat(sub.hScores);
XX = repmat(x, size(YY,1), 1);
RR = (rand(size(XX))-0.5)/10;
% [~,idx] = sort(ss);
% imagesc(XX(idx,:))

axes(ax)
% figure(1); clf;
hold on;
% plot(x,vertcat(sub.hScores)', '.','color', [1 1 1]*0.75)
l1 = plot(RR'+XX',YY', '-','color', [1 1 1]*0.75);
ff = 0.5;
l2 = plot(RR'+XX',YY', 'o','color', AZblue*ff+(1-ff), 'markersize', 5);
l3 = plot(x,nanmean(vertcat(sub.hScores)), '-','color', AZblue, 'linewidth', 6);
xlim([0.75 4.25])

set(gca, 'xtick', [1:4], 'tickdir', 'out', 'fontsize', 18)
xlabel('suspicion score')
ylabel('frequency')