function weights = epan(bw)
    weights = (-bw:bw) ./ bw;
    weights = arrayfun(@(k) 1-k^2, weights);
    sumweights = sum(weights);
    weights = (weights ./ sumweights)';
end


