%% ECO 513 - Problem Set 4
% Ryan Broll

%% Q1

% (ii) Degree of invertibility
    
    % Compute inputs to Kalman Filter
        % H, F, Q, and E_00
n = 6;
H = [1 0 -.3 2 .8 0]';
F = [0 0 0 0 0 0; 0 0 0 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0];
Q = zeros(6);
Q(1,1) = 1;
Q(2,2) = 1;
E_00 = zeros(6,1);

        % P_00
vecP_00 = (eye(n^2) - kron(F,F))\(Q(:));
I_n = eye(n);
P_00 = kron(I_n(:)',I_n)*kron(I_n,vecP_00);

        % Generate y from DGP to filter
T = 10000;
eps = randn(T+2,2);
y = zeros(T,1);
for t = 1:T
y(t) = eps(t+2,1) - .3*eps(t+1,1) + 2*eps(t+1,2) + .8*eps(t,1);
end

    % Compute E_tt and P_tt using Kalman filter
[E_tt, P_tt, E_tT, P_tT] = kalman_filter(H,F,Q,E_00,P_00,y);

    % Display degrees of invertibility
disp('*****Q1 (ii)******')
disp('eps_1 Degree of invertibility')
disp(1 - P_tt(1,1))
disp('eps_2 Degree of invertibility')
disp(1 - P_tt(2,2))

% (iii) Degree of recoverability
    % Display degrees of recoverability
disp('*****Q1 (iii) *****')
disp('eps_1 Degree of recoverability')
disp(1 - P_tT(1,1))
disp('eps_2 Degree of recoverability')
disp(1 - P_tT(2,2))


%% Q3

% (i) SVAR-IV IRF
    
    % Fix parameters: p = # of lags, hor = # of horizons, n = dimension of
    % y
p = 4;
hor = 36;
n = 4;

    % Load data
q3table = readtable('gk_data.csv');
y = q3table{:,["gs1" "ebp" "ipgr" "infl"]};

    % Estimate VAR(p) on y
[Ahat, Sigmahat, etahat] = var_estim(y,p);

    % Calculate Hcol1, the first column of H
z = q3table{p+1:end,"ff4_tc"};
gammahat = mean(etahat .* z)';
Hcol1 = (1/sqrt((gammahat'/Sigmahat)*gammahat))*gammahat;

    % Recursively compute reduced form impulse responses Psihat
Ahat_3d = reshape(Ahat(:,2:end),n,n,p); % 3D version of Ahat (just slopes)
Psihat = zeros(n,n,hor+1); % Also a 3D array
Psihat(:,:,1) = eye(n); % Psihat_0 = I_n

for l = 1:p % Calculate first p reduced form impulse responses
    Psihat(:,:,l+1) = sum(pagemtimes(Ahat_3d(:,:,1:l),Psihat(:,:,l:-1:1)),3);
end

for l = p+1:hor % Calculate the rest
    Psihat(:,:,l+1) = sum(pagemtimes(Ahat_3d, Psihat(:,:,l:-1:l-3)),3);
end

    % Compute actual impulse responses
Thetahat = zeros(hor+1,n);
for l = 0:hor
    Thetahat(l+1,:) = (Psihat(:,:,l+1) * Hcol1)';
end

    % Normalize Thetahat
Thetahat_normalized = Thetahat/Thetahat(1,1);

    % Plot
plot(0:hor, Thetahat_normalized(:,2))
xlabel("horizon (months)")
ylabel("percent")
title("Impulse Response of EBP to Monetary Shock: SVAR implementation")

% (ii) SVAR-IV FVD

fvd = zeros(hor,1);
Sigmahatrep = repmat(Sigmahat,1,1,hor+1);
for l=1:hor
    num = Thetahat(1:l,2)'*Thetahat(1:l,2);
    denom = pagemtimes(pagemtimes(Psihat,Sigmahatrep),'none',Psihat,'transpose');
    denom = sum(denom(:,:,1:l),3);
    fvd(l) = num/denom(1,1);
end

    %Plot
plot(1:hor, fvd)
xlabel('horizon l (months)')
ylabel('FVD_{ebp,1,l}')
title("Forecast Variance Decomposition of EBP with respect to Monetary Shock")


%% Q3 continued

% (iii) LP-IV: recursive VAR implementation

    % Estimate VAR(4) with z ordered first
z = q3table{:,"ff4_tc"};
[Ahat3iii, Sigmahat3iii, uhat3iii, X3iii] = var_estim([z y], p);

    % Recursively compute reduced form impulse responses Psihat
Ahat3d_3iii = reshape(Ahat3iii(:,2:end),n+1,n+1,p); % 3D version of Ahat (just slopes)
Psihat3iii = zeros(n+1,n+1,hor+1); % Also a 3D array
Psihat3iii(:,:,1) = eye(n+1); % Psihat_0 = I_n

for l = 1:p % Calculate first p reduced form impulse responses
    Psihat3iii(:,:,l+1) = sum(pagemtimes(Ahat3d_3iii(:,:,1:l),Psihat3iii(:,:,l:-1:1)),3);
end

for l = p+1:hor % Calculate the rest
    Psihat3iii(:,:,l+1) = sum(pagemtimes(Ahat3d_3iii, Psihat3iii(:,:,l:-1:l-3)),3);
end

    % Calculate actual impulse responses using recursive identification
H = chol(Sigmahat3iii, 'lower');

Theta3iii = zeros(hor+1,n+1);
for l = 1:hor+1
    ir3iii = Psihat3iii(:,:,l) * H;
    Theta3iii(l,:) = ir3iii(:,1)';
end

    % Normalize
Theta3iii_normalized = Theta3iii / Theta3iii(1,2);

    % Plot
plot(0:hor, Theta3iii_normalized(:,3))
xlabel("horizon (months)")
ylabel("percent")
title("Impulse Response of EBP to Monetary Shock: LP-IV, VAR implementation")

% (v) Test of invertibility - assuming i.i.d. Gaussian shocks
p = 4;
n = 5;
k = n*p+1;
r = 16;
Omegahat3iii = kron(inv(X3iii'*X3iii), Sigmahat3iii);
R=eye(n*k);
R = R([7:10 32:35 57:60 82:85],:);
Wald_stat = (R*Ahat3iii(:))' * ( (R*Omegahat3iii*R') \ (R*Ahat3iii(:)) );
p_value = 1-chi2cdf(Wald_stat,r);

%%
% (iv) LP-IV: 2SLS implementation
z = q3table{:,"ff4_tc"};
y = q3table{:,["gs1" "ebp" "ipgr" "infl"]};

beta_l = zeros(hor+1,1);
for l = 0:hor
    effT = length(z) - p - l;
    y_iv = y(p+1+l:end,2);
    lagcontrols = lagmatrix([y z],1:p);
    X_iv = [ones(effT,1) y(p+1:end-l,1) lagcontrols(p+1:end-l,:)];
    Z_iv = [ones(effT,1) z(p+1:end-l) lagcontrols(p+1:end-l,:)];
    W_iv = inv(Z_iv'*Z_iv);
    S_ZX = Z_iv'*X_iv;
    S_ZY = Z_iv'*y_iv;
    denom = S_ZX'*W_iv*S_ZX;
    betahat = denom\(S_ZX'*W_iv*S_ZY); % GMM estimator
    beta_l(l+1) = betahat(2);
end

    % Plot
plot(0:hor, beta_l)
xlabel("horizon (months)")
ylabel("percent")
title("Impulse Response of EBP to Monetary Shock: LP-IV, 2SLS implementation")


%% Question 5

% Load Data

    % INDPRO
indpro = table2timetable(readtable('INDPRO.csv'));
indpro = synchronize(indpro, lag(indpro));
indpro.indpro = 4*log(indpro.INDPRO_indpro ./ indpro.INDPRO_2);
indpro = indpro(:,'indpro');

    % PCEPI
pcepi = table2timetable(readtable('PCEPI.csv'));
pcepi = synchronize(pcepi, lag(pcepi));
pcepi.pcepi = 4*log(pcepi.PCEPI_pcepi ./ pcepi.PCEPI_2);
pcepi = pcepi(:,'pcepi');

    % GS1
gs1 = table2timetable(readtable('GS1.csv'));
gs1 = synchronize(gs1, pcepi);
gs1.gs1 = .01*gs1.GS1 - gs1.pcepi;
gs1 = gs1(:,'gs1');

    % Combine and impose sample restrictions
table = synchronize(indpro,pcepi,gs1);
sample = timerange('1984-01-01', '2020-01-01');
table = table(sample,:);
data = table{:,:};

% (i) Bounds for identified set
p = 12;
n= 3;
[Ahat5, Sigmahat5] = var_estim(data,p);

    % Recursively compute reduced form impulse responses Psihat
Ahat5_3d = reshape(Ahat5(:,2:end),n,n,p); % 3D version of Ahat (just slopes)
Psihat5 = zeros(n,n,hor+1); % Also a 3D array
Psihat5(:,:,1) = eye(n); % Psihat_0 = I_n

for l = 1:p % Calculate first p reduced form impulse responses
    Psihat5(:,:,l+1) = sum(pagemtimes(Ahat5_3d(:,:,1:l),Psihat5(:,:,l:-1:1)),3);
end

for l = p+1:hor % Calculate the rest
    Psihat5(:,:,l+1) = sum(pagemtimes(Ahat5_3d, Psihat5(:,:,l:-1:l-(p-1))),3);
end

    % Generate distribution of impulse responses 
C = chol(Sigmahat5)';
numdraws = 1000000;

impulse_dist = zeros(hor+1,2); %[min .16 .5 .84 max] for each l = 0,...,36

for l = 0:hor
disp(l)
admissQs = sign_restrict(Psihat5,C,n,numdraws);
PsiC_rep = repmat(Psihat5(1,:,l+1)*C,1,1,size(admissQs,3));
impulse = pagemtimes(PsiC_rep, admissQs(:,1,:));
impulse_dist(l+1,1) = min(impulse);
impulse_dist(l+1,5) = max(impulse);
impulse_dist(l+1,2:4) = quantile(impulse,[.16 .5 .84]);
end

%%
plot(0:36,impulse_dist(:,1), "r")
hold on
plot(0:36,impulse_dist(:,2), "g")
plot(0:36,impulse_dist(:,3), "k")
plot(0:36,impulse_dist(:,4), "g")
plot(0:36,impulse_dist(:,5), "r")
title("Output Growth Impulse Response: Identified and 68% Credible Set")
legend('ID set','68% band', 'estimate')
xlabel('horizon (months)')
ylabel('impulse response')
hold off




