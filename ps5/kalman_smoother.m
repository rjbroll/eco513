function [S_smooth, P_smooth] = kalman_smoother(F, S_filter, P_filter, S_forec, P_forec)

    % Kalman smoother
    % Trans. eqn:   S_t = F*S_{t-1} + epsilon_t
    
    
    %% Recursion

    T = size(S_filter,1);
    S_smooth = zeros(size(S_filter));
    P_smooth = zeros(size(P_filter));
    
    % Initialize using filter at t=T
    S_smooth(end,:) = S_filter(end,:);
    P_smooth(:,:,end) = P_filter(:,:,end);
    
    for t=T-1:-1:1
        
        J_t = P_filter(:,:,t)*(F'/P_forec(:,:,t+1));
        S_smooth(t,:) = S_filter(t,:) + (S_smooth(t+1,:)-S_forec(t+1,:))*J_t';
        P_smooth(:,:,t) = P_filter(:,:,t) + J_t*(P_smooth(:,:,t+1)-P_forec(:,:,t+1))*J_t';
        
    end

end