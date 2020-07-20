function [phitScore, phitN] = scorePHIT_v1(phit, phitEmails)

EM = vertcat(phit.emails);
CL = vertcat(phit.click);

for i = 1:length(phitEmails)
    ind = strcmp(EM, phitEmails{i});
    phitScore(i,1) = nanmean(CL(ind));
    phitN(i,1) = sum(ind);
end
