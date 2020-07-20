function [sitEmails, ty, id] = get_sitEmails_v1(sub)

% get unique set of sitEmails used in a group of subjects
% ty - label for type (+1 real vs 0 simulated)
% id - identify (+1 phish vs 0 safe)

EM = vertcat(sub.emailCode);
TY = strcmp(vertcat(sub.type), 'Type'); % 1 for real
ID = strcmp(vertcat(sub.realId), 'Scam'); % 1 for scam
[sitEmails, i1, i2] = unique(EM);
ty = TY(i1); id = ID(i1);