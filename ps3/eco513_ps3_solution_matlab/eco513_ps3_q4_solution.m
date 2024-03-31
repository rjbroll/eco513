clear;

% ECO 513, Spring 2024
% Solutions to PS3, Q4
% Mikkel Plagborg-Moller, 2024-03-18


%% Settings

% VAR specification
p = 4; % VAR lags

% Bootstrap
numrep = 1e4; % Number of bootstrap replications
rng(20240318); % Set random number seed


%% Load data

dat_monthly = table2timetable(readtable('eco513_ps1_Monthly.txt'));
dat_quarterly = table2timetable(readtable('eco513_ps1_Quarterly.txt'));
dat = innerjoin(retime(dat_monthly,'quarterly','mean'), dat_quarterly);

rgdp_growth = 4*log(dat{2:end,'GDPC1'}./dat{1:end-1,'GDPC1'}); % RGDP log growth, annualized
pgdp_growth = 4*log(dat{2:end,'GDPDEF'}./dat{1:end-1,'GDPDEF'}); % PGDP log growth, annualized
fedfunds_real = 0.01*dat{2:end,'FEDFUNDS'}-pgdp_growth; % Fed Funds Rate minus inflation

% Keep relevant sample
time = datetime(dat.Properties.RowTimes(2:end));
sample = (time >= datetime('1954-07-01')) & (time <= datetime('2019-10-01'));
time = time(sample);
rgdp_growth = rgdp_growth(sample);
pgdp_growth = pgdp_growth(sample);
fedfunds_real = fedfunds_real(sample);

Y = [rgdp_growth pgdp_growth fedfunds_real]; % VAR variables;
[T,n] = size(Y); % VAR dimensions


%% Estimate VAR

inds_test = 3*(1:p)'; % Indices of coefficients corresponding to lags of PGDP growth
[wald_stat, betahat, Sigmahat, res] = var_waldstat(Y, p, 1, inds_test); % Wald statistic and VAR residuals


%% Residual bootstrap

disp('Bootstrapping...')
tic;

% Construct parameters that impose null hypothesis
betahat_boot = betahat;
betahat_boot(1,inds_test) = 0;

wald_stats_boot = zeros(numrep, 1);

for i=1:numrep % For each bootstrap replication...
   
    % Draw p initial observations at random from data
    the_init_ind = ceil(rand*(T-p+1));
    the_Y_init = Y(the_init_ind:the_init_ind+p-1,:);
    
    % Bootstrap residuals
    the_res = res(randi(T-p,T-p,1),:);
    
    % Recursively generate bootstrap data
    the_Y = [the_Y_init; zeros(T-p,n)];
    
    for t=p+1:T
        the_Y(t,:) = betahat_boot(:,1)' + the_res(t-p,:); % Constant plus bootstrap residuals
        for l=1:p % Add lagged terms
            the_Y(t,:) = the_Y(t,:) + the_Y(t-l,:)*betahat_boot(:,2+(l-1)*n:1+l*n)';
        end
    end
    
    % Compute bootstrap Wald statistic
    wald_stats_boot(i) = var_waldstat(the_Y, p, 1, inds_test);
    
    % Print progress
    if mod(i, floor(numrep/10))==0
        fprintf('%3d%s\n', 100*i/numrep, '%');
    end
    
end


%% Report results

disp('Elapsed time (seconds)');
disp(toc);

disp('Wald statistic');
disp(wald_stat);

disp('Bootstrap quantiles of null distribution [0.5 0.75 0.9 0.95 0.99]');
disp(quantile(wald_stats_boot, [0.5 0.75 0.9 0.95 0.99]));

disp('Bootstrap p-value');
disp(mean(wald_stats_boot > wald_stat));


