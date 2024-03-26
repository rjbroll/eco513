function [E_tt, P_tt, E_tT, P_tT] = kalman_filter(H,F,Q,E_00,P_00,y)

% Initialize filter vectors
T = length(y);
dim = size(H,1);
E = zeros(dim,T+1);
P = zeros(dim,dim,T+1);
E_next = zeros(dim,T);
P_next = zeros(dim,dim,T);
E(:,1) = E_00;
P(:,:,1) = P_00;

% Run Kalman filter
for i = 1:T
    E_next(:,i) = F*E(:,i);
    P_next(:,:,i) = F*P(:,:,i)*F' + Q;
    y_next = H'*E_next(:,i);
    h_t = H'*P_next(:,:,i)*H;
    K_t = P_next(:,:,i)*(H/(h_t));
    eta_t = y(i) - y_next;
    E(:,i+1) = E_next(:,i) + K_t*eta_t;
    P(:,:,i+1) = P_next(:,:,i) - K_t*H'*P_next(:,:,i);
end

E_tt = E(:,T+1);
P_tt = P(:,:,T+1);

% Initialize smoother vectors
E_smooth = zeros(dim,T);
P_smooth = zeros(dim,dim,T);
E_smooth(:,T) = E_tt;
P_smooth(:,:,T) = P_tt;

% Run Kalman smoother
for i = T-1:1
    J_t = P(:,:,i+1)*(F'/(P_next(:,:,i+1)));
    E_smooth(:,i) = E(:,i+1) + J_t*(E_smooth(:,i+1) - E_next(:,i+1));
    P_smooth(:,:,i) = P(:,:,i+1) + J_t*(P_smooth(:,:,i+1) - P_next(:,:,i+1))*J_t';
end

E_tT = E_smooth(:,1);
P_tT = P_smooth(:,:,1);

end