function B = run_regressionModel_v1(sub, nPast)

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
