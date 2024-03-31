function restricted_Qs = sign_restrict(Psi, C, n, numdraws)
    % Draw from Haar measure
    Q = haar_draws(n,numdraws);

    % Check sign restrictions
    pass_check = zeros(numdraws,1); % 1 if Q passes sign check, 0 if not
    for i = 1:numdraws
        check = zeros(n,1);
        for l = 0:2
        check(l+1) = ((Psi(3,:,l+1)*C*Q(:,1,i) >= 0) & (Psi(2,:,l+1)*C*Q(:,1,i) <= 0));
        end
        pass_check(i) = (sum(check) == n);
    end

    % Extract Qs that pass the check
    restricted_Qs = Q(:,:,logical(pass_check));
end