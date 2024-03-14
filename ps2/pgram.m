function [omegagrid, pgram] = pgram(y)

    % set up omega grid
    T = length(y);
    jgrid = -floor((T-1)/2):floor(T/2);
    omegagrid = ((2*pi)/T).*jgrid;

    % compute fft
    ffty = fft(y)';

    % adjust fft to create periodogram
    zindex = floor((T+1)/2); % index of omega=0
    evec = arrayfun(@(k) exp(-1j*k), omegagrid);
    part1 = (1/sqrt(T)).*ffty(T-zindex+2:T).*evec(1:zindex-1);
    part2 = (1/sqrt(T)).*ffty(1:T-zindex+1).*evec(zindex:T);
    pgram = arrayfun(@(k) abs(k)^2, [part1 part2]);
    pgram(zindex) = pgram(zindex+1); % fix bias at 0

end