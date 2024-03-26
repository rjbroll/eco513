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
title("Impulse Response of Excess Bond Premium to Monetary Shock")

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




















