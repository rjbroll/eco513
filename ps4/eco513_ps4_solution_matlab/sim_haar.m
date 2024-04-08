function Q = sim_haar(n, numsim)

    % Simulate from Haar measure

    Q = zeros(n,n,numsim);
    for i=1:numsim
        X = randn(n);
        [the_Q,the_R] = qr(X);
        Q(:,:,i) = the_Q.*sign(diag(the_R)');
    end

end