rng(1)
clear
m = 30;
nSigs = 6;
masterInput = randn(m, nSigs);
dp = randn(nSigs, 1);
nSections = 3;
polesOrig = rand(m * 5 * nSections, 1);
%% My implementation
maspoles = [polesOrig; masterInput(:)];
for idx = 1 : nSections
    initialValueZ2 = randn(1, nSigs) * 0.1;
    initialValueZ1 = randn(1, nSigs) * 0.1;
    maspoles = [maspoles; initialValueZ2(:); initialValueZ1(:);];
end
poles = ADNode(maspoles);
flat = poles(m * 5 * (nSections - 1) + m * 5 + 1 : m * 5 * (nSections - 1) + m * 5 + m * nSigs);
input = reshape(flat, size(masterInput));
coeff = cell(nSections, 1);
states = cell(nSections, 1);
for idx = 1 : nSections
    base = m * 5 * (idx - 1);
    coeff{idx}.b0 = poles(base + 1 : base + m);
    coeff{idx}.d1 = poles(base + m + 1 : base + m * 2);
    coeff{idx}.d2 = poles(base + m * 2 + 1 : base + m * 3);
    coeff{idx}.c1 = poles(base + m * 3 + 1 : base + m * 4);
    coeff{idx}.c2 = poles(base + m * 4 + 1 : base + m * 5);
    base = m * 5 * (nSections - 1) + m * 5 + m * nSigs + nSigs * 2 * (idx - 1);
    states{idx}.z2_A = poles(base + 1 : base + nSigs);
    states{idx}.z2_A = reshape(states{idx}.z2_A, fliplr(size(states{idx}.z2_A.value)));
    states{idx}.z1_A = poles(base + nSigs + 1 : base + nSigs * 2);
    states{idx}.z1_A = reshape(states{idx}.z1_A, fliplr(size(states{idx}.z1_A.value)));
end
% Forward filtering
for a = 1 : m - 1
    out = input(a, :);
    for idx = 1 : nSections
        y1 = out - states{idx}.z1_A - states{idx}.z2_A;
        out = coeff{idx}.b0(a) .* y1 + coeff{idx}.d1(a) .* states{idx}.z1_A + coeff{idx}.d2(a) .* states{idx}.z2_A;
        states{idx}.z2_A = states{idx}.z2_A + coeff{idx}.c2(a) .* states{idx}.z1_A;
        states{idx}.z1_A = states{idx}.z1_A + coeff{idx}.c1(a) .* y1;
    end
    if a == 1
        tmp2 = out;
    else
        tmp2(a, :) = out;
    end
end
out = input(m, :);
for idx = 1 : nSections
    y1 = out - states{idx}.z1_A - states{idx}.z2_A;
    out = coeff{idx}.b0(m) .* y1 + coeff{idx}.d1(m) .* states{idx}.z1_A + coeff{idx}.d2(m) .* states{idx}.z2_A;
end
tmp2(m, :) = out;
sub = tmp2 * dp;
mse1 = sum(sub.^2);
grad1 = mse1.backprop(1);
input_grad = grad1(m * 5 * (nSections - 1) + m * 5 + 1 : m * 5 * (nSections - 1) + m * 5 + m * nSigs);
input_grad = reshape(input_grad, size(masterInput));
coeff_grad = cell(nSections, 1);
states_grad = cell(nSections, 1);
for idx = 1 : nSections
    base = m * 5 * (idx - 1);
    coeff_grad{idx}.b0 = grad1(base + 1 : base + m);
    coeff_grad{idx}.d1 = grad1(base + m + 1 : base + m * 2);
    coeff_grad{idx}.d2 = grad1(base + m * 2 + 1 : base + m * 3);
    coeff_grad{idx}.c1 = grad1(base + m * 3 + 1 : base + m * 4);
    coeff_grad{idx}.c2 = grad1(base + m * 4 + 1 : base + m * 5);
    base = m * 5 * (nSections - 1) + m * 5 + m * nSigs + nSigs * 2 * (idx - 1);
    states_grad{idx}.z2_A = grad1(base + 1 : base + nSigs);
    states_grad{idx}.z2_A = reshape(states_grad{idx}.z2_A, fliplr(size(states_grad{idx}.z2_A)));
    states_grad{idx}.z1_A = grad1(base + nSigs + 1 : base + nSigs * 2);
    states_grad{idx}.z1_A = reshape(states_grad{idx}.z1_A, fliplr(size(states_grad{idx}.z1_A)));
end
plot(grad1);axis tight