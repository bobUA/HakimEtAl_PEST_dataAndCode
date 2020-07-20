function get_demographics_v1(sub)

% number of women
disp(['n women = ' num2str(sum(strcmp({sub.sex}', 'F')))])
disp(['n men = ' num2str(sum(strcmp({sub.sex}', 'M')))])

% mean age
disp(['mean age = ' num2str(nanmean([sub.age]))]);
disp(['max age = ' num2str(max([sub.age]))]);
disp(['min age = ' num2str(min([sub.age]))]);
disp(['std age = ' num2str(nanstd([sub.age]))]);
