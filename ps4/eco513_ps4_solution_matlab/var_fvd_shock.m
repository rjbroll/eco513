function fvds = var_fvd_shock(beta, Sigma, Hcol, maxhorz)

    % VAR forecast variance decomposition
    
    n = size(beta,1);
    
    % Impulse responses
    [irfs, irfs_rf] = var_irf_shock(beta, Hcol, maxhorz);
    
    % Calculate FVD
    fvds = zeros(n,maxhorz);
    numer = zeros(n,1);
    denom = zeros(n,1);
    for l=1:maxhorz
        numer = numer + irfs(:,l).^2;
        denom = denom + diag(irfs_rf(:,:,l)*Sigma*irfs_rf(:,:,l)');
        fvds(:,l) = numer./denom;
    end

end