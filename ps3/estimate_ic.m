function [ic_matrix, aic_selectedp, bic_selectedp]  = estimate_ic(y, p_range)

% set size parameters
T = size(y,1) ;
n = size(y,2) ;
length_p_range = length(p_range) ;
p_max = max(p_range) ;
effT = T - p_max ;

% pre-allocate IC matrix. First column is p, 2nd is aic, 3rd is bic.
ic_matrix = zeros(length_p_range,3) ;

% estimate aic, bic for each p in p_range
for i = 1:length_p_range
    p = p_range(i) ;
    ic_matrix(i,1) = p ;
    
    [~,Sigmahat] = estimate_var(y((p_max+1 - p):end,:),p) ;
    ldSh = log(det(Sigmahat)) ;
    ic_matrix(i,2) = ldSh + n^2*p*2/effT ;
    ic_matrix(i,3) = ldSh + n^2*p*log(effT)/effT ;
end

[~, aic_bestindex] = max(-ic_matrix(:,2)) ;
[~, bic_bestindex] = max(-ic_matrix(:,3)) ;
aic_selectedp = p_range(aic_bestindex) ;
bic_selectedp = p_range(bic_bestindex) ;

end
