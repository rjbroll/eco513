function [p_aic, p_bic, aics, bics] = var_ic(Y, p_min, p_max)

    % AIC/BIC for VAR(p)
    
    % Compute log(det(Sigmahat)) for different models
    logdets = zeros(p_max-p_min+1,1);
    for p=p_min:p_max
        [~, the_Sigmahat] = var_estim(Y(p_max-p+1:end,:),p); % Enforce same estimation sample size
        logdets(p-p_min+1) = log(det(the_Sigmahat));
    end
    
    % Compute AIC/BIC criteria
    T = size(Y,1)-p_max;
    n = size(Y,2);
    
    aics = logdets + n^2*(p_min:p_max)'*2/T;
    bics = logdets + n^2*(p_min:p_max)'*log(T)/T;
    
    % Find minimum
    [~,minind_aic] = min(aics);
    [~,minind_bic] = min(bics);
    
    p_aic = minind_aic+p_min-1;
    p_bic = minind_bic+p_min-1;

end