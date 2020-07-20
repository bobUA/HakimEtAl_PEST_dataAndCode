function [score, lab] = computeScores_v1(sub)

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
