function p = dGamma(x, r, lambda)

p = lambda^r * x.^(r-1) .* exp(-lambda * x) ./ gamma(r);