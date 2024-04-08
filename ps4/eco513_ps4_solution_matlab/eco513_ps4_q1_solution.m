clear;

% ECO 513, Spring 2024
% Solutions to PS4, Q1
% Mikkel Plagborg-Moller, 2024-04-04


%% Settings

% Lag lengths to consider in Kalman filter/smoother ...
% when computing degree of invertibility and recoverability
data_lengths = [100 200 300];


%% Model

% MA coefficients in the form
% Y_t = sum_{j=0}^{n_eps} sum_{l=0}^q Theta_{j,l} epsilon_{j,t-l}
Theta = [1  -0.3   0.8;
             0  2      0];

% Set up Kalman matrices
% State vector: S_t = (epsilon_{1,t}, ... epsilon_{n_eps,t}, epsilon_{1,t-1}, ..., epsilon_{n_eps,t-q})'
[n_eps,qp1] = size(Theta);
c = 0; % Intercept
H = Theta(:)'; % Measurement matrix
F = [zeros(n_eps,n_eps*qp1);
     eye(n_eps*(qp1-1)), zeros(n_eps*(qp1-1),n_eps)]; % Transition matrix
Q = blkdiag(eye(n_eps), zeros(n_eps*(qp1-1))); % Transition var-cov
R = 0; % Measurement var-cov
S_1_0 = zeros(n_eps*qp1,1); % Initial mean
P_1_0 = eye(n_eps*qp1); % Initial variance


%% Degree of invertibility/recoverability

for T=data_lengths 
   
    disp('Data vector length');
    disp(T);
    
    % Run filter
    [~, S_filter, P_filter, S_forec, P_forec] = kalman_filter(zeros(T,1), c, H, F, Q, R, S_1_0, P_1_0); % Feed in zero data
    disp('Degree of invertibility for each shock');
    disp(1-diag(P_filter(1:n_eps,1:n_eps,T))');
    
    % Run smoother
    [~, P_smooth] = kalman_smoother(F, S_filter, P_filter, S_forec, P_forec);
    disp('Degree of recoverability for each shock');
    disp(1-diag(P_smooth(1:n_eps,1:n_eps,floor(T/2)))');
    
    disp(' ');
    
end


