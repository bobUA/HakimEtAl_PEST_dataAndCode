function [sitScore, sitN, sitRT] = scoreSIT_v1(sub, refEmails)

% score SIT based on list of emails
% NOTE: has some hacks to deal with different coding of emails between PHIT
% and SIT data sets

emails = vertcat(sub.emailCode);
rating = vertcat(sub.userId);
RT = vertcat(sub.RT);
for i = 1:length(refEmails)
    ind = strcmp(emails, refEmails{i});
    sitScore(i,1) = nanmean(rating(ind));
    sitRT(i,1) = nanmedian(RT(ind));
    sitN(i,1) = sum(ind);
end
