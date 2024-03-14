function [omegagrid, spec_est, spec_low, spec_high] = spec_ar(y, p)

% Estimate VAR with p_adj lags
[betahat, sigmahat, uhat, X] = estimate_var(y,p);
ahat = betahat(:,2:(p+1))'; % just slope coefs - no intercept

% Set up grid of omega values
effT = length(y) - p;
jgrid = 0:floor(effT/2);
omegagrid = ((2*pi)/effT).*jgrid;

% Compute estimate of the spectrum
spec_est = arrayfun(@(omega) spec_at_omega(omega, ahat, sigmahat), omegagrid)

% Compute pointwise 95% CI bounds for spectrum estimate using delta method
seval = arrayfun(@(omega) se_at_omega(omega, ahat, sigmahat, uhat, X), omegagrid)
spec_low = spec_est - 1.96.*seval;
spec_high = spec_est + 1.96.*seval;
end


% function that computes log spectrum at specific omega grid value
function specval = spec_at_omega(omega, ahat, sigmahat)
evec = arrayfun(@(l) exp(-1j*omega*l),1:length(ahat));
specval = log(sigmahat) - log(2*pi) - log(abs(1- evec*ahat)^2);
end

% function that computes standard error at specific omega grid value
function seval = se_at_omega(omega, ahat, sigmahat, uhat, X)

% set constants
effT = length(uhat);
p = length(ahat);

% calculate A matrix - 'meat' of delta-method sandwich
V_ut2 = ((uhat.^2 - sigmahat)'*(uhat.^2 - sigmahat))/(effT); % Var(u_t^2)
Sx_inv = inv(X'*X)./effT;
Omegahat = Sx_inv .* sigmahat;
Omegahat_minus1 = Omegahat(2:end, 2:end); % asymp. cov. matrix of just slope coefs
zero = zeros(1,p);
A = [V_ut2 zero; zero' Omegahat_minus1];

% calculate B matrix - 'bread' of delta-method sandwich
evec = arrayfun(@(l) exp(-1j*omega*l),1:length(ahat));
dvec = arrayfun(@(q) dfunc(q, omega, ahat), 1:p);  % 'inside' derivatives wrt each slope coef
outside = abs(1-(evec*ahat))^2;  % 'outside' derivative
B = [1/sigmahat dvec./outside];

% calculate standard error for spectrum at omega value
seval = sqrt(B*A*B'/effT);
end


% function that computes inside derivative wrt qth slope coef
function dq = dfunc(q, omega, ahat);
p = length(ahat);
cos_vec = arrayfun(@(l) cos(omega*(q-l)), 1:p);
dq = -2*cos(omega*q) + 2*cos_vec*ahat;
end



