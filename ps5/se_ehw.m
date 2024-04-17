function stderr_ehw = se_ehw(X, eps, coefinds)
    Xeps = eps .* X;
    V_ehw = (X'*X)\Xeps' * Xeps/(X'*X);
    stderr_ehw = diag(V_ehw);
    stderr_ehw = sqrt(stderr_ehw(coefinds)); 
end