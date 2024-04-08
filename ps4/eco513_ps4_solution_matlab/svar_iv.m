function [irfs, fvds, betahat, res] = svar_iv(Y, Z, p, maxhorz)

    % SVAR-IV point estimation (single IV)
    
    [betahat, Sigmahat, res] = var_estim(Y, p);
    gammahat = res'*Z(p+1:end)/length(res);
    Hcolhat = gammahat/sqrt(gammahat'*(Sigmahat\gammahat));

    % IRFs
    irfs = var_irf_shock(betahat, Hcolhat, maxhorz);

    % FVDs
    fvds = var_fvd_shock(betahat, Sigmahat, Hcolhat, maxhorz);

end