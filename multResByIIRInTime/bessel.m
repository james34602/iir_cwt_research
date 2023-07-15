function [b, a] = bessel(n, timeConstant)
nP1 = n + 1;
i = 0 : n;
%% Compute pole
a = zeros(1, nP1);
if timeConstant <= eps
    a(1) = 1;
    b = 1;
    return;
end
for idx = 1 : nP1
    k = idx - 1;
    a(idx) = (-1)^k * binomial(n, k) * prod((2 * timeConstant + i) ./ (2 * timeConstant + i + k));
end
%% Compute unity gain compensation
y = filter(1, [1 -1], a);
b = y(nP1);
end
function nk = binomial(n, k)
if k == 0
    nk = 1;
    return;
else
    nk = n;
end
for i = 2 : k
    nk = (nk * (n - (i - 1))) / i;
end
end