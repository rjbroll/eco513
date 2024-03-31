clear;

% ECO 513, Spring 2024
% Solutions to PS3, Q2
% Mikkel Plagborg-Moller, 2024-03-18


%% Settings

% Estimation settings
Z_lags = 3; % Number of IV lags of inflation and unemployment rate
p_max = 20; % Max number of lags for BIC in VAR-HAC routine
hac_fct = @(Y) varhac(Y, p_max); % HAC function


%% Load data

load_data; % Reads from data files and enforces common sample
data_matrix_lag = lagmatrix(data_matrix_sample,1:Z_lags);

% Response variable: inflation
Y = data_matrix_sample(Z_lags+1:end-1,1);

% Regressors: constant, future inflation, lagged inflation, unemployment rate
T = length(Y);
X = [ones(T,1) data_matrix_sample(Z_lags+2:end,1) data_matrix_sample(Z_lags:end-2,1) data_matrix_sample(Z_lags+1:end-1,2)];

% Instruments: constant, lags of inflation and unemployment rate
Z = [ones(T,1) data_matrix_lag(Z_lags+1:end-1,:)];

% Dimensions
k = size(X,2);
r = size(Z,2);


%% GMM estimation

% 1st step
[betahat_1st, betahat_var_1st, Omegahat_1st] = gmm_iv(Y, X, Z, inv(Z'*Z), hac_fct, false);

% 2nd (efficient) step
[betahat_eff, betahat_var_eff, ~, J_stat] = gmm_iv(Y, X, Z, inv(Omegahat_1st), hac_fct, true);

% J test p-value
J_pval = 1-chi2cdf(J_stat, r-k);


%% Report results

disp('1st step: point estimates [const pi(t+1) pi(t-1) x(t)]');
disp(betahat_1st');

disp('1st step: standard errors');
disp(sqrt(diag(betahat_var_1st)'));

disp('2nd step: point estimates');
disp(betahat_eff');

disp('2nd step: standard errors');
disp(sqrt(diag(betahat_var_eff)'));

disp('J test: statistic');
disp(J_stat);

disp('J test: p-value');
disp(J_pval);

