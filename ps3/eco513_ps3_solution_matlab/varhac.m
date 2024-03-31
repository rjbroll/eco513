function [Omegahat, p_bic] = varhac(Y, p_max)

    % VAR-HAC estimator, using BIC to select lag length

    [~,p_bic] = var_ic(Y,0,p_max); % Obtain p by BIC
    [betahat, Sigmahat] = var_estim(Y,p_bic); % Estimate VAR(p)
    
    n = size(Y,2);
    Ahat_reshape = reshape(betahat(:,2:end),n,n,p_bic); % Reshape hat{A} into 3-dim array
    
    % VAR(p) long-run variance
    denom = eye(n)-sum(Ahat_reshape,3);
    Omegahat = (denom\Sigmahat)/(denom');

end