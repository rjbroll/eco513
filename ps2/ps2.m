%% Q2 

% Simulate 10,000 AR time-series with T = 200
simdata = draws(10000, 200, 0, .6, .2, 1);

% Estimate p_hat by AIC and BIC, estimate rho_1 using p_hat, repeat 10,000
% times

% pre-allocate results matrix. First column is aic, second is bic.
phat_matrix = zeros(10000,2) ;
rhohat_matrix = zeros(10000,2) ;

for r = 1:10000
    % Estimate IC
    [~, aic_p, bic_p] = estimate_ic(simdata(r,:)',1:3) ;
    phat_matrix(r,:) = [aic_p, bic_p] ;

    % Estimate VAR(selected p) for AIC
    aic_betahat = estimate_var(simdata(r,(4-aic_p):end)', aic_p) ;
    rhohat_matrix(r,1) = aic_betahat(2) ;

    % Estimate VAR(selected p) for BIC
    bic_betahat = estimate_var(simdata(r,(4-bic_p):end)',bic_p) ;
    rhohat_matrix(r,2) = bic_betahat(2) ;
end 

%% Graph rho_1_hat

    % Graph of rho_1_hat sampling density
[aicrhoy, aicrhox] = ksdensity(rhohat_matrix(:,1)) ;
[bicrhoy, bicrhox] = ksdensity(rhohat_matrix(:,2)) ;
figure; 
plot(aicrhox, aicrhoy, bicrhox, bicrhoy) ;
legend('AIC', 'BIC') ;
title('Sampling Density of rho_1 hat') ;

%% Graph p_hat

    % Graph of p_hat sampling density
figure
histogram(phat_matrix(:,1), Normalization ='pdf')
xticks(1:3)
ylim([0,1])
yticks(0:.2:1)
title('AIC')


figure
histogram(phat_matrix(:,2), Normalization = 'pdf')
xticks(1:3)
ylim([0,1])
yticks(0:.2:1)
title('BIC')


%% Q5

% Import data
un_adj = readtable("UNRATE.csv");
un_not = readtable("UNRATENSA.csv");
data_adj = [un_adj.UNRATE];
data_not = [un_not.UNRATENSA];

T = length(data_adj);


%%% (i)(a) Estimate log spectral density for adjusted data by AR(p) %%%

% Estimate p for data_adj by BIC
[~,~,p_adj] = estimate_ic(data_adj,1:50);

% Estimate spectral density and pointwise 95% CI bounds
[omega_grid_var_adj, spec_est_var_adj, spec_low_var_adj, spec_high_var_adj] = spec_ar(data_adj, p_adj);

%%
figure
plot(omega_grid_var_adj, spec_high_var_adj, 'r')
hold on
plot(omega_grid_var_adj, spec_est_var_adj, 'k')
hold on
plot(omega_grid_var_adj, spec_low_var_adj, 'r')
xlabel('frequency')
ylabel('log spectrum')
title('Seasonally adjusted')





%%
%%% (i)(b) Estimate log spectral density for unadjusted data by AR(p) %%%

% Estimate p for data_not by BIC
[~,~,p_not] = estimate_ic(data_not, 1:50);

% Estimate spectral density and pointwise 95% CI bounds
[omega_grid_var_not, spec_est_var_not, spec_low_var_not, spec_high_var_not] = spec_ar(data_not, p_not);

%%
figure
plot(omega_grid_var_not, spec_high_var_not, 'r')
hold on
plot(omega_grid_var_not, spec_est_var_not, 'k')
hold on
plot(omega_grid_var_not, spec_low_var_not, 'r')
xlabel('frequency')
ylabel('log spectrum')
title('Not seasonally adjusted')


%%
%%% (ii)(a) Estimate log spectral density for adjusted data by kernel %%%

% Compute the periodogram
[omegagrid_2a, pgram_adj] = pgram(data_adj);
log_pgram_adj = log(pgram_adj);

% smooth the periodogram
T = length(log_pgram_adj);
weights = epan(10);
wide_pgram_adj = [log_pgram_adj log_pgram_adj log_pgram_adj]; % create wide version to smooth boundaries
smooth_pgram_adj = zeros(1,T);
for i = 1:T
    smooth_pgram_adj(i) = wide_pgram_adj(T+i-10 : T+i+10) * weights;
end
smooth_pgram_adj = smooth_pgram_adj(floor((T+1)/2):T);

% compute 95% CI bounds
se = sqrt(sum(weights.^2));
smooth_pgram_adj_high = smooth_pgram_adj + 1.96*se;
smooth_pgram_adj_low = smooth_pgram_adj - 1.96*se;
xgrid = omegagrid_2a(floor((T+1)/2):T);

%% plot
figure
plot(xgrid, smooth_pgram_adj, 'k')
hold on
plot(xgrid, smooth_pgram_adj_low, 'r')
hold on
plot(xgrid, smooth_pgram_adj_high, 'r')
xlabel('frequency')
ylabel('log spectrum')
title('Seasonally adjusted')


%%
%%% (ii)(b) Estimate log spectral density for unadjusted data by kernel %%%

% Compute the periodogram
[omegagrid_2b, pgram_not] = pgram(data_not);
log_pgram_not = log(pgram_not);

% smooth the periodogram
T = length(log_pgram_not);
weights = epan(10);
wide_pgram_not = [log_pgram_not log_pgram_not log_pgram_not]; % create wide version to smooth boundaries
smooth_pgram_not = zeros(1,T);
for i = 1:T
    smooth_pgram_not(i) = wide_pgram_not(T+i-10 : T+i+10) * weights;
end
smooth_pgram_not = smooth_pgram_not(floor((T+1)/2):T);

% compute 95% CI bounds
se = sqrt(sum(weights.^2));
smooth_pgram_not_high = smooth_pgram_not + 1.96*se;
smooth_pgram_not_low = smooth_pgram_not - 1.96*se;
xgrid = omegagrid_2b(floor((T+1)/2):T);

%% plot
figure
plot(xgrid, smooth_pgram_not, 'k')
hold on
plot(xgrid, smooth_pgram_not_low, 'r')
hold on
plot(xgrid, smooth_pgram_not_high, 'r')
xlabel('frequency')
ylabel('log spectrum')
title('Not seasonally adjusted')







