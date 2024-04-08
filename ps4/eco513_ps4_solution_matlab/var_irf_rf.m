function irfs = var_irf_rf(beta, maxhorz)

    % Reduced-form impulse responses for estimated VAR

    n = size(beta,1);
    A = reshape(beta(:,2:end),n,n,[]); % VAR coefficients
    p = size(A,3);
    
    % Recursion to compute IRFs
    irfs = zeros(n,n,maxhorz+1);
    irfs(:,:,1) = eye(n);
    for l=1:maxhorz
        for j=1:min(l,p)
            irfs(:,:,l+1) = irfs(:,:,l+1) + A(:,:,j)*irfs(:,:,l-j+1);
        end
    end

end