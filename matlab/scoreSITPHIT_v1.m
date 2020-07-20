function [sitScore, sitN] = scoreSITPHIT_v1(sub, phitEmails)

% score SIT based on PHIT emails
% NOTE: has some hacks to deal with different coding of emails between PHIT
% and SIT data sets

emails = vertcat(sub.emailCode);
rating = vertcat(sub.userId);
for i = 1:length(phitEmails)
    if length(phitEmails{i})  == 5
        em = ['0' phitEmails{i}];
    elseif length(phitEmails{i})  == 4
        em = ['00' phitEmails{i}];
    else
        em = phitEmails{i};
    end
    ind = strcmp(emails, em);
    sitScore(i,1) = nanmean(rating(ind));
    sitN(i,1) = sum(ind);
end
