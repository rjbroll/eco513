function Y_sim = var_boot(beta, Y, res)

    % Recursive residual VAR bootstrap

    % Dimensions
    [n,npp1] = size(beta);
    p = (npp1-1)/n;
    T = size(res,1);
    
    % Bootstrap residuals
    res_boot = res(ceil(rand(T-p,1)*(T-p)),:);
    
    % Draw initial value
    init_ind = ceil(rand()*(T-p));
    Y_init = Y(init_ind:init_ind+p,:);
    
    % Simulate recursively
    Y_sim = [Y_init; zeros(T-p, n)];
    for t=p+1:T
        Y_sim(t,:) = beta(:,1)' + res_boot(t-p,:);
        for l=1:p
            Y_sim(t,:) = Y_sim(t,:) + Y_sim(t-l,:)*beta(:,2+(l-1)*n:1+l*n)';
        end
    end

end