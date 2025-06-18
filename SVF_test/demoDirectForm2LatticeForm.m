addpath('../')
[b, a] = butter(1, [0.1, 0.2], 'stop');
b0=b(1);b1=b(2);b2=b(3);a2=a(2);a3=a(3);
%% TF2LATTIC
[k, v] = tf2latc(b, a);
k = [0; k];
kDelta = [1; zeros(128, 1)];
exampleOut = filter(b, a, kDelta);
N = length(v);
g = zeros(N, 1);
y = zeros(length(kDelta), 1);
for n = 1 : length(kDelta)
    f = kDelta(n);
    acc = 0;
    for m = N : -1 : 2
        f = f - k(m) * g(m - 1);
        g(m) = k(m) * f + g(m - 1);
        acc = acc + v(m) * g(m);
    end
    g(1) = f;
    y(n) = acc + v(1) * g(1);
end