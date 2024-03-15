function [betahat, Sigmahat, uhat, X] = estimate_var(y, p)   

    % create regressor matrix X (effT x k)
    X = lagmatrix(y,1:p) ;
    effT = size(y,1) - p ;
    X = [ones(effT,1) X(p+1:end,:)] ;

    % create y (effT x n)
    y = y(p+1:end,:) ;

    % estimate betahat (nxk)
    betahat = (X\y)' ;

    % estimate Sigmahat (nxn)
    uhat = y - X*betahat' ;
    Sigmahat = (uhat'*uhat)/effT ;


end