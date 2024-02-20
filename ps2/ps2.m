%% Q2 %%%%%%%%%%%%%%%%

% Simulate 10,000 AR time-series with T = 200
simdata = draws(10000, 200, 0, .6, .2, 1);

% Estimate p_hat by AIC and BIC, estimate rho_1 using p_hat, repeat 10,000
% times

% pre-allocate results matrix. First column is aic, second is bic.
phat_matrix = zeros(10000,2) ;
rhohat_matrix = zeros(10000,2) ;

for r = 1:10000
    % Estimate IC
    [~, aic_p, bic_p] = estimate_ic(simdata(r,:)',1:3) ;
    phat_matrix(r,:) = [aic_p, bic_p] ;

    % Estimate VAR(selected p) for AIC
    aic_betahat = estimate_var(simdata(r,(4-aic_p):end)', aic_p) ;
    rhohat_matrix(r,1) = aic_betahat(2) ;

    % Estimate VAR(selected p) for BIC
    bic_betahat = estimate_var(simdata(r,(4-bic_p):end)',bic_p) ;
    rhohat_matrix(r,2) = bic_betahat(2) ;
end 

%%
% Graph results
    % Graph of rho_1_hat sampling density
[aicrhoy, aicrhox] = ksdensity(rhohat_matrix(:,1)) ;
[bicrhoy, bicrhox] = ksdensity(rhohat_matrix(:,2)) ;
figure; 
plot(aicrhox, aicrhoy, bicrhox, bicrhoy) ;
legend('AIC', 'BIC') ;
title('Sampling Density of rho_1 hat') ;

%%

    % Graph of p_hat sampling density
figure
histogram(phat_matrix(:,1), Normalization ='pdf')
xticks(1:3)
ylim([0,1])
yticks(0:.2:1)
title('AIC')


figure
histogram(phat_matrix(:,2), Normalization = 'pdf')
xticks(1:3)
ylim([0,1])
yticks(0:.2:1)
title('BIC')












%% Q5 %%%%%%%%%%%%%%%%%%%%%%%%





