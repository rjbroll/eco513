function [irfs, betahat, res] = svar_chol(Y, p, maxhorz)

    % Cholesky-identified SVAR

    [betahat, Sigmahat, res] = var_estim(Y, p); % Estimate augmented VAR
    Hhat = chol(Sigmahat, 'lower'); % Cholesky decomposition of H
    irfs = var_irf_shock(betahat, Hhat(:,1), maxhorz); % IRFs to first shock

end
