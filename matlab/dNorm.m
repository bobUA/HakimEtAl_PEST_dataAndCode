function p = dNorm(x, mu, tau)

% note tau is precision
p = sqrt(tau/2/pi) * exp(-tau * (x - mu).^2 / 2);