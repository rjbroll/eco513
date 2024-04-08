clear;

% ECO 513, Spring 2024
% Solutions to PS4, Q5
% Mikkel Plagborg-Moller, 2024-04-04


%% Settings

p = 12;                 % Lag length for VAR
K = 2;                  % Max horizon at which to impose sign restrictions
maxhorz = 36;           % Max horizon to compute identified set for
numsim = 1e6;           % Number of simulated Q matrices
alpha = 0.32;           % Significance level for "posterior" credible band


%% Load data

dat = readtable('eco513_ps4_Monthly.txt');
ipgrowth = 400*log(dat.INDPRO(2:end)./dat.INDPRO(1:end-1)); % Industrial production growth
infl = 400*log(dat.PCEPI(2:end)./dat.PCEPI(1:end-1)); % Inflation
rir = dat.GS1(2:end) - infl; % Real interest rate

Y = [rir infl ipgrowth]; % Data matrix

[T,n] = size(Y);


%% Part (i)

% Estimate reduced-form VAR and compute relevant parameters
[betahat, Sigmahat] = var_estim(Y, p);
C = chol(Sigmahat)'; % Sigmahat=C*C'
irfs_rf = var_irf_rf(betahat, maxhorz); % Reduced-form IRFs

% Set up linear inequalities M*Q>=0 that impose sign restrictions
M = zeros(2*K+2,3);
signs = [1 -1];
icol = 0;
for i=1:2
    for l=0:K
        icol = icol+1;
        M(icol,:) = signs(i)*irfs_rf(i,:,l+1)*C;
    end
end

% Simulate Qs that satisfy the sign restrictions
Qs = sim_haar(n, numsim); % Draw several Qs from Haar measure
Q1s = permute(Qs(:,1,:), [1 3 2]); % Keep only first column
Q1s_satisfy = Q1s(:, all(M*Q1s>=0)); % Keep only those satisfying sign restrictions

% Compute IRF of industrial production wrt. these Qs
irf_draws = permute(irfs_rf(3,:,:), [3 2 1])*C*Q1s_satisfy;


%% Part (ii)

% Quantiles of IRF draws
quants = quantile(irf_draws,[alpha/2 0.5 1-alpha/2],2);


%% Part (iii)

figure;
plot(0:maxhorz, [min(irf_draws,[],2) max(irf_draws,[],2)], '-k'); % Lines: estimate of identified set
hold on;
fill([0:maxhorz fliplr(0:maxhorz)], [quants(:,1)' fliplr(quants(:,3)')], [0.7 0.7 0.7], 'EdgeColor', 'none'); % Shaded area between quantiles
plot(0:maxhorz, quants(:,2), '-r'); % "Posterior" median
hold off;
xlabel('horizon (months)');
title('output response (annualized percent)');
xlim([0 maxhorz]);
set(gca, 'XTick', 0:6:maxhorz);

