function theta = min_dist(vecA)

    % Function mapping VAR coefficients vec(A) into structural parameters theta

    rho_pi = vecA(1);
    rho_x = vecA(3);

    gamma_f = 1 / (rho_pi + 1 / (1 + vecA(2)/vecA(4)*rho_x/(1-rho_pi)));
    lambda = (1 - gamma_f*rho_pi)*rho_x / vecA(4) - gamma_f*rho_x;

    theta = [gamma_f; lambda; rho_pi; rho_x];

end