clear;

% ECO 513, Spring 2024
% Solutions to PS3, Q3
% Mikkel Plagborg-Moller, 2024-03-18


%% Settings

% Finite difference settings
stepsize = 1e-10; % Step size


%% Load data

load_data; % Reads from data files and enforces common sample

% VAR(1) estimation
[betahat, Sigmahat, res, X] = var_estim(data_matrix_sample, 1);
betahat_var = kron(inv(X'*X), Sigmahat);
Ahat = betahat(:,2:end); % VAR(1) coefficient matrix hat{A}
Psihat = betahat_var(3:end,3:end); % Var-cov matrix of hat{A}


%% Compute minimum distance estimate

thetahat = min_dist(Ahat(:)); % hat{theta}


%% Compute standard errors

% Obtain derivatives by finite differences
k = length(thetahat);
Rhatinv = zeros(k);
the_I = eye(k);
for j=1:k
    thetahat_plus = min_dist(Ahat(:) + stepsize*the_I(:,j));
    thetahat_minus = min_dist(Ahat(:) - stepsize*the_I(:,j));
    Rhatinv(:,j) = (thetahat_plus - thetahat_minus)/(2*stepsize); % dtheta/dmu_j, mu = vec(A)
end

thetahat_var = Rhatinv*Psihat*Rhatinv'; % Var-cov matrix of hat{theta}


%% Report results

disp('Estimates [gamma_f lambda rho_pi rho_x]');
disp(thetahat');

disp('Standard errors');
disp(sqrt(diag(thetahat_var)'));




