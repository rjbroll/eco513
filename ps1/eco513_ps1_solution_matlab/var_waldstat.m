function [wald_stat, betahat, Sigmahat, res] = var_waldstat(Y, p, eqn, inds_test)
    
    % Compute Wald statistic for a test that certain coefficients are zero
    % in a particular equation of a VAR(p) model
    % assuming homoskedastic innovations

    [betahat, Sigmahat, res, X] = var_estim(Y, p); % VAR(p) estimation
    betahat_eqn_var = inv(X'*X)*Sigmahat(eqn,eqn); % Variance of coefficient estimates in given equation

    % Wald statistic
    wald_stat = betahat(eqn,inds_test)*(betahat_eqn_var(inds_test,inds_test)\betahat(eqn,inds_test)');

end