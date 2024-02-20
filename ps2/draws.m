
% function to simulate R AR time series given parameters
function simdata = var_sim(R, T, nu, rho1, rho2, sigma2)
    simdata = zeros(R,T+2) ;
    u_matrix = normrnd(0,sigma2,R,T) ;
    for r = 1:R
        for t = 3:(T+2)
            simdata(r,t) = nu + rho1*simdata(r,t-1) + rho2*simdata(r,t-2) + u_matrix(r,t-2) ;
        end
    end
    simdata = simdata(:,3:end) ;
end
