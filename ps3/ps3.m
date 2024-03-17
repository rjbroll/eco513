

%% Q2

% Load data
    % UNRATE - take quarterly averages
unrate = table2timetable(readtable('UNRATE.csv'));
unrate = retime(unrate, 'quarterly','mean');

    % GDPDEF - take annualized log growth rates
deflator = table2timetable(readtable('GDPDEF.csv'));
deflator = synchronize(deflator, lag(deflator));
deflator.growth = 4*log(deflator.GDPDEF_deflator ./ deflator.GDPDEF_2);
deflator = deflator(:,'growth');

    % Combine
table = innerjoin(deflator, unrate);

    % Select sample
sample = timerange('-inf', '2020-01-01');
table = table(sample,:);

    % Turn into data matrix (instead of timetable)
data = table{:,["growth" "UNRATE"]};

% (i) Estimate betahat by GMM
    % Set effective sample size
effT = size(data,1) - 3 - 1; % extra -1 because of the x_{t+1}

    % Create effT x 7 matrix of instruments Z
Z_1 = lagmatrix(data(:,1),1:3);
Z_2 = lagmatrix(data(:,2),1:3);
Z = [ones(effT,1) Z_1(4:end-1,:) Z_2(4:end-1,:)]; 

    % Create effT x 4 matrix of "regressors" X
X = lagmatrix(data(:,1),[-1 1]);
X = [ones(effT,1) X(4:end-1,:) data(4:end-1,2)];

    % Create effT x 1 vector y
y = data(4:end-1,1);

    % Estimate betahat with GMM (W = inv(Z'Z) equivalent to 2SLS)
fitX = Z*(Z\X); % fitted values from regression of X on Z
fity = Z*(Z\y); % fitted values from regression of y on Z
betahat = fitX\fity; % GMM = 2SLS 

% (ii) Two-step efficient GMM
    % Create effT x 7 matrix M of GMM moment observations
M = (y-X*betahat).*Z;

    % Estimate the LRV Omegahat of the time series process M by VAR-HAC
        % Choose p using BIC
[~,~,bic_p] = estimate_ic(M,1:50);

        % Estimate VAR(bic_p) on M
[A,Sigmahat] = estimate_var(M,bic_p);

        % Construct Omegahat
A = A(:,2:end); % keep only "slope" coefficient estimates 
I = repmat(eye(7),bic_p,1);
sumA = A * I;
Omegahat = ((eye(7)-sumA)\Sigmahat)/((eye(7)-sumA)');

    % Estimate betahat with GMM efficient (W = inv(Omegahat))
betahat_eff = (((X'*Z)/Omegahat)*(Z'*X))\(((X'*Z)/Omegahat)*(Z'*y));

% (iii) Standard errors for (i) and (ii)
    % Estimate 7x4 matrix Ghat, the derivative of moments wrt parameters
Ghat = -(Z'*X)/effT;

    % Standard errors for (i): Avar1 = inv(G'WG)*G'WOmegaWG*inv(G'WG),
    % W=inv(Z'Z)
outer = ((Ghat')/(Z'*Z))*Ghat;
inner = ((Ghat')/(Z'*Z)) * Omegahat * ((Z'*Z)\Ghat);
AVar1 = (outer\(inner))/outer;
se1 = sqrt(diag(AVar1)/effT);

    % Standard errors for (ii): Avar2 = inv(G'inv(Omega)G)
AVar2 = inv((Ghat'/Omegahat)*Ghat);
se2 = sqrt(diag(AVar2)/effT);

% (iv) Hansen's J-test
    % Create matrix M_eff of GMM moment observations using betahat_eff
M_eff = (y-X*betahat_eff).*Z;
    % Construct test statistic and run test
Qhat_eff = (mean(M_eff,1)/(Omegahat)) * mean(M_eff,1)';
Jstat = effT * Qhat_eff;
pvalue = 1 - chi2cdf(Jstat,3); % dof = r-k = 7-4 = 3


%% Q3

% (ii) Estimate theta by minimum distance
    % run a VAR(1) on (x_t, pi_t)'
[Ahat_q3, Sigmahat_q3, ~, X_q3] = estimate_var(fliplr(data),1);

    % Construct Psihat, the asymptotic variance of Ahat_q3
effT_q3 = size(data,1) - 1;
Sx = (X_q3'*X_q3)/effT_q3;
Psihat = kron(inv(Sx), Sigmahat_q3);
Psihat = Psihat(3:end,3:end); % just variance of slope coefficients

    % Conduct CMD with efficient weight matrix W = inv(Psihat)
Ahat_q3 = (Ahat_q3(:,2:end));
Ahat_q3 = Ahat_q3(:); % Ahat_q3 = [A11 A12 A21 A22]'
objfun = @(theta) (((Ahat_q3 - h(theta))')/Psihat) * (Ahat_q3 - h(theta));
[thetahat, opval] = fminunc(objfun, [.56,-.01,.76,-.01]');

% (iii) Calculate standard errors
    % Calculate D
D = zeros(4);

theta_1p = thetahat + [abs(thetahat(1))/100 0 0 0]';
theta_1m = thetahat - [abs(thetahat(1))/100 0 0 0]';
D(:,1) = (h(theta_1p) - h(theta_1m))/(abs(thetahat(1))/50);

theta_2p = thetahat + [0 abs(thetahat(2))/100 0 0]';
theta_2m = thetahat - [0 abs(thetahat(2))/100 0 0]';
D(:,2) = (h(theta_2p) - h(theta_2m))/(abs(thetahat(2))/50);

theta_3p = thetahat + [0 0 abs(thetahat(3))/100 0]';
theta_3m = thetahat - [0 0 abs(thetahat(3))/100 0]';
D(:,3) = (h(theta_3p) - h(theta_3m))/(abs(thetahat(3))/50);

theta_4p = thetahat + [0 0 0 abs(thetahat(4))/100]';
theta_4m = thetahat - [0 0 0 abs(thetahat(4))/100]';
D(:,4) = (h(theta_4p) - h(theta_4m))/(abs(thetahat(4))/50);

    % Calculate standard errors
AVar_q3 = inv((D'/Psihat)*D);
se_q3 = sqrt(diag(AVar_q3)/effT_q3);














