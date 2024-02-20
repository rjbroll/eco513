function [betahat, Sigmahat] = estimate_var(y, p)   

    % create regressor matrix
    X = lagmatrix(y,1:p) ;
    effT = size(y,1) - p ;
    X = [ones(effT,1) X(p+1:end,:)] ;

    % create y
    y = y(p+1:end,:) ;

    % estimate betahat
    betahat = (X\y)' ;

    % estimate Sigmahat
    uhat = y - X*betahat' ;
    Sigmahat = (uhat'*uhat)/effT ;

end