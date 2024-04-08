function betahat = ivreg_2sls(Y, X, Z, R)

    % 2SLS regression of Y on X, using IVs Z and controls R
    
    % Data matrices
    [n,k] = size(X);
    Rt = [ones(n,1) R]; % Add intercept
    Xt = [X Rt];
    Zt = [Z Rt];
    
    % Estimate
    aux = (Xt'*Zt)/(Zt'*Zt);
    aux2 = (aux*(Zt'*Xt))\aux;
    betahat_grand = aux2*(Zt'*Y);
    betahat = betahat_grand(1:k);

end