function irf = lpiv_2sls(Y, Z, p, resp_ind, norm_ind, maxhorz)

    % LP-IV IRF estimation, 2SLS implementation
    
    irf = zeros(1,1+maxhorz);
    W_lag = lagmatrix([Y Z],1:p);
    
    for l=0:maxhorz
        irf(l+1) = ivreg_2sls(Y(p+l+1:end,resp_ind), ...    % Outcome
                                Y(p+1:end-l,norm_ind), ...  % Endogenous regressor
                                Z(p+1:end-l), ...           % IV
                                W_lag(p+1:end-l,:));        % Controls
    end

end