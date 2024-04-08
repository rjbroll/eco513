clear;

% ECO 513, Spring 2024
% Solutions to PS4, Q3
% Mikkel Plagborg-Moller, 2024-04-04


%% Settings

p = 4;                  % Lag length for VARs and LPs
maxhorz = 36;           % Maximum IRF/FVD horizon

num_boot = 1e3;         % Number of bootstrap draws
signif = 0.1;           % Significance level

plot_ind = 4;           % Plot IRFs and FVDs for this variable
plot_ylim_irf = [-0.5 2]; % y-axis limits for IRF plots


%% Load data

dat = readtable('gk_data.csv');
Y_names = {'gs1', 'ipgr', 'infl', 'ebp'};
Y = dat{:,Y_names}; % First variable will be normalized to have impact response 1 to shock
Z = dat{:,'ff4_tc'}; % IV
W = [Z Y]; % Joint data matrix
[T,n] = size(Y);


%% Parts (i)-(ii)

% Point estimates
[irfs_svariv, fvds_svariv, betahat_svariv, res_svariv] = svar_iv(Y, Z, p, maxhorz); % SVAR-IV
irf_svariv_norm = irfs_svariv(plot_ind,:)/irfs_svariv(1,1); % Normalize scale and sign of shock

% Construct matrix beta = [c, A_1, ..., A_p] for expanded VAR in ...
% (Z_t, Y_t')' that treats Z_t as i.i.d.
betahat_svariv_wZ = zeros(n+1,1+p*(n+1));
betahat_svariv_wZ(:,1) = [mean(Z); betahat_svariv(:,1)]; % Intercept
for l=1:p
    betahat_svariv_wZ(:,2+(l-1)*(n+1):1+l*(n+1)) = [zeros(1,n+1); zeros(n,1) betahat_svariv(:,2+(l-1)*n:1+l*n)]; % Lag matrices (zeros for IV)
end

% Run bootstrap
irf_svariv_norm_boot = zeros(num_boot, 1+maxhorz);
fvd_svariv_boot = zeros(num_boot, maxhorz);
disp('Bootstrapping SVAR-IV');

for b=1:num_boot
    
    the_W_sim = var_boot(betahat_svariv_wZ, W, [Z(p+1:end) res_svariv]); % Bootstrap data
    [the_irfs_svariv, the_fvds_svariv] = svar_iv(the_W_sim(:,2:end), the_W_sim(:,1), p, maxhorz); % Estimate
    
    % Store results
    irf_svariv_norm_boot(b,:) = the_irfs_svariv(plot_ind,:)/the_irfs_svariv(1,1);
    fvd_svariv_boot(b,:) = the_fvds_svariv(plot_ind,:);
    
    % Print progress
    if mod(b, ceil(num_boot/10))==0
        fprintf('%3d%s\n', round(100*b/num_boot), '%');
    end
    
end

% Plot IRF
q_boot_irf_svariv_norm = quantile(irf_svariv_norm_boot, [signif/2,1-signif/2]); % Bootstrap quantiles
plot_band(0:maxhorz, irf_svariv_norm, ...
          q_boot_irf_svariv_norm(1,:), q_boot_irf_svariv_norm(2,:), ...
          'SVAR-IV: IRF', plot_ylim_irf);

% Plot FVD
q_boot_fvds_svariv = quantile(fvd_svariv_boot, [signif/2,1-signif/2]); % Bootstrap quantiles
plot_band(1:maxhorz, fvds_svariv(plot_ind,:), ...
          q_boot_fvds_svariv(1,:), q_boot_fvds_svariv(2,:), ...
          'SVAR-IV: FVD', [0 1]);


%% Parts (iii)-(iv)

% VAR estimates
[irfs_lpiv_var, betahat_lpiv_var, res_lpiv_var] = svar_chol(W, p, maxhorz); % Cholesky SVAR
irf_lpiv_var_norm = irfs_lpiv_var(1+plot_ind,:)/irfs_lpiv_var(2,1); % Normalize

% 2SLS point estimates
irf_lpiv_2sls = lpiv_2sls(Y, Z, p, plot_ind, 1, maxhorz);

% Run bootstrap
irf_lpiv_var_norm_boot = zeros(num_boot, 1+maxhorz);
irf_lpiv_2sls_boot = zeros(num_boot, 1+maxhorz);
disp('Bootstrapping LP-IV');

for b=1:num_boot
    
    the_W_sim = var_boot(betahat_lpiv_var, W, res_lpiv_var); % Bootstrap data
    the_irfs_lpiv_var = svar_chol(the_W_sim, p, maxhorz); % VAR estimate
    the_irfs_lpiv_2sls = lpiv_2sls(the_W_sim(:,2:end), the_W_sim(:,1), p, plot_ind, 1, maxhorz); % 2SLS estimate
    
    % Store results
    irf_lpiv_var_norm_boot(b,:) = the_irfs_lpiv_var(1+plot_ind,:)/the_irfs_lpiv_var(2,1);
    irf_lpiv_2sls_boot(b,:) = the_irfs_lpiv_2sls;
    
    % Print progress
    if mod(b, ceil(num_boot/10))==0
        fprintf('%3d%s\n', round(100*b/num_boot), '%');
    end
    
end

% Plot VAR estimate
q_boot_irf_lpiv_var_norm = quantile(irf_lpiv_var_norm_boot, [signif/2,1-signif/2]); % Bootstrap quantiles
plot_band(0:maxhorz, irf_lpiv_var_norm, ...
          q_boot_irf_lpiv_var_norm(1,:), q_boot_irf_lpiv_var_norm(2,:), ...
          'LP-IV, VAR implementation: IRF', plot_ylim_irf);

% Plot 2SLS estimate
q_boot_irf_lpiv_2sls = quantile(irf_lpiv_2sls_boot, [signif/2,1-signif/2]); % Bootstrap quantiles
plot_band(0:maxhorz, irf_lpiv_2sls, ...
          q_boot_irf_lpiv_2sls(1,:), q_boot_irf_lpiv_2sls(2,:), ...
          'LP-IV, 2SLS implementation: IRF', plot_ylim_irf);


%% Part (v)

disp('Granger causality p-values');
for i=1:n
    disp(Y_names{i});
    the_wald_stat = var_waldstat(W, p, i+1, 2+(0:p-1)*(n+1)); % Test that coefs on IV are zero in equation i+1
    disp(1-chi2cdf(the_wald_stat,p)); % p-value
end
