function [ll, S_filter, P_filter, S_forec, P_forec] = kalman_filter(Y, c, H, F, Q, R, S_1_0, P_1_0)

    % Kalman filter
    % Obs. eqn:     Y_t = c + H*S_t + U_t
    % Trans. eqn:   S_t = F*S_{t-1} + epsilon_t
    % Shocks:       U_t ~ N(0,R), epsilon_t ~ N(0,Q)
    
    
    %% Initialize
    
    S_t_tm1 = S_1_0;
    P_t_tm1 = P_1_0;
    [T,n] = size(Y);
    m = length(S_1_0);
    
    
    %% Recursion
    
    ll = 0;
    S_filter = zeros(T,m);
    P_filter = zeros(m,m,T);
    S_forec = zeros(T,m);
    P_forec = zeros(m,m,T);
    
    for t=1:T
        
        if nargout>1
            S_forec(t,:) = S_t_tm1';
            P_forec(:,:,t) = P_t_tm1;
        end
       
        % Forecast Y
        Y_t_tm1 = c + H*S_t_tm1;
        V_t_tm1 = H*P_t_tm1*H' + R;
        forec_err = Y(t,:)' - Y_t_tm1;
        
        % Log likelihood contribution
        ll = ll - 0.5*(n*log(2*pi) + log(det(V_t_tm1)) + forec_err'*(V_t_tm1\forec_err));
        
        % Update S
        M_t = P_t_tm1*(H'/V_t_tm1);
        S_t_t = S_t_tm1 + M_t*forec_err;
        P_t_t = (eye(m) - M_t*H)*P_t_tm1;
        if nargout>1
            S_filter(t,:) = S_t_t';
            P_filter(:,:,t) = P_t_t;
        end
        
        % Forecast S
        S_t_tm1 = F*S_t_t;
        P_t_tm1 = F*P_t_t*F' + Q;
        
    end

end