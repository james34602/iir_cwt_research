addpath('../')
nSigmas = 6;
sigmas = [100, 200, 150, 20, 40, 3];
[b, a, c1, c2] = gauss_precompute(sigmas);
coeff = cell(length(sigmas), 1);
for idx = 1 : length(sigmas)
    [k, v] = tf2latc(b(idx), a(idx, :));
    N = length(v);
    coeff{idx}.k = [0; k];
    coeff{idx}.v = v;
end
g = zeros(N, 1);
%%
sel = 1;
dcSignal = ones(10000, 1);
oo1 = zeros(size(dcSignal, 1), 1);
oo2 = zeros(size(dcSignal, 1), 1);
oo3 = zeros(size(dcSignal, 1), 1);
z1_A = 0.0;
z2_A = 0.0;
z0 = 0; z1 = 0;
s1 = 0; s2 = 0;
for idx = 1 : size(dcSignal, 1)
    if idx >= 1500
        if mod(idx, 500) == 0
            sel = sel + 1;
            if sel > length(sigmas)
                sel = 1;
            end
        end
    end
    %% Coeff update
    x = dcSignal(idx);
    %% Simplified SVF for diagram
    %     y = x - z1_A - z2_A;
    %     oo1(idx) = b0(sel) * y + 2 * c2(sel) * z1_A + z2_A;
    %     z2_A = z2_A + c2(sel) * z1_A;
    %     z1_A = z1_A + c1(sel) * y;
    %% Simplified SVF for runtime
    y = x - z1_A - z2_A;
    st = c2(sel) * z1_A;
    oo1(idx) = b(sel) * y + 2 * st + z2_A;
    z2_A = z2_A + st;
    z1_A = z1_A + c1(sel) * y;
    %% Direct form
    Yi = b(sel) * x + z0;
    z0 = z1 - a(sel, 2) * Yi;
    z1 = -a(sel, 3) * Yi;
    oo2(idx) = Yi;
    %% Lattice form
    f = x;
    acc = 0;
    for m = N : -1 : 2
        f = f - coeff{sel}.k(m) * g(m - 1);
        g(m) = coeff{sel}.k(m) * f + g(m - 1);
        acc = acc + coeff{sel}.v(m) * g(m);
    end
    g(1) = f;
    oo3(idx) = acc + coeff{sel}.v(1) * g(1);
end
figure(1)
plot(oo1);
hold on
plot(oo2);
plot(oo3);
hold off
axis tight;
legend('State variable', 'Direct form', 'Lattice')