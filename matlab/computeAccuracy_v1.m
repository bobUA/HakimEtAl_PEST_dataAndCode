function [ACC, AA, safeScore, scamScore] = computeAccuracy_v1(sub)


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