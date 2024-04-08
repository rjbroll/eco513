function [irfs, irfs_rf] = var_irf_shock(beta, Hcol, maxhorz)
    
    % Rotation of reduced-form VAR impulse responses
    
    n = size(beta,1);
    
    % Reduced-form IRFs
    irfs_rf = var_irf_rf(beta, maxhorz);
    
    % Reshape and rotate
    irfs = zeros(n, maxhorz+1);
    for l=0:maxhorz
        irfs(:,l+1) = irfs_rf(:,:,l+1)*Hcol;
    end

end