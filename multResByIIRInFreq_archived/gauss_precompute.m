function [b, a, c1, c2] = gauss_precompute(sigma)
sigma = sigma(:);
fnc1 = @(x) ((606 * x.^2) / 1087 - (3009 * x) / 5513 + 712 / 5411) ./ (x.^2 - (496 * x) / 541 + 719 / 1034);
mp = exp(-137 ./ (100 * sigma));
a1 = 2 * cos(fnc1(sigma) ./ sigma) .* mp;
%% Transfer function
a = [ones(size(sigma, 1), 1), -a1, mp .* mp];
b = sum(a, 2);
%% SVF realization
c1 = 2 - a1;
c2 = b ./ c1;
end