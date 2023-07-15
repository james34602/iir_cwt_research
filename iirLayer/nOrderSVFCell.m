rng(1)
clear
m = 256;
nSigs = 16;
masterInput = randn(m, nSigs);
dp = randn(nSigs, 1);
nSections = 4;
polesOrig = rand(m * 5 * nSections, 1);
target = randn(m, 1);
reverseFilter = 1;
%% My implementation
poles = [polesOrig; masterInput(:)];
for idx = 1 : nSections
    initialValueZ2 = randn(1, nSigs) * 0.1;
    initialValueZ1 = randn(1, nSigs) * 0.1;
    poles = [poles; initialValueZ2(:); initialValueZ1(:);];
end
flat = poles(m * 5 * (nSections - 1) + m * 5 + 1 : m * 5 * (nSections - 1) + m * 5 + m * nSigs);
input = reshape(flat, size(masterInput));
coeff = cell(nSections, 1);
states = cell(nSections, 1);
misc = cell(nSections, 1);
for idx = 1 : nSections
    base = m * 5 * (idx - 1);
    coeff{idx}.b0 = poles(base + 1 : base + m);
    coeff{idx}.d1 = poles(base + m + 1 : base + m * 2);
    coeff{idx}.d2 = poles(base + m * 2 + 1 : base + m * 3);
    coeff{idx}.c1 = poles(base + m * 3 + 1 : base + m * 4);
    coeff{idx}.c2 = poles(base + m * 4 + 1 : base + m * 5);
    base = m * 5 * (nSections - 1) + m * 5 + m * nSigs + nSigs * 2 * (idx - 1);
    states{idx}.z2_A = poles(base + 1 : base + nSigs);
    states{idx}.z2_A = reshape(states{idx}.z2_A, fliplr(size(states{idx}.z2_A)));
    states{idx}.z1_A = poles(base + nSigs + 1 : base + nSigs * 2);
    states{idx}.z1_A = reshape(states{idx}.z1_A, fliplr(size(states{idx}.z1_A)));
    misc{idx}.y2 = zeros(1, nSigs);
    misc{idx}.z2_A = states{idx}.z2_A;
    misc{idx}.z1_A = states{idx}.z1_A;
end
if reverseFilter
    input = flipud(input);
    for idx = 1 : nSections
        coeff{idx}.b0 = flipud(coeff{idx}.b0);
        coeff{idx}.d1 = flipud(coeff{idx}.d1);
        coeff{idx}.d2 = flipud(coeff{idx}.d2);
        coeff{idx}.c1 = flipud(coeff{idx}.c1);
        coeff{idx}.c2 = flipud(coeff{idx}.c2);
    end
end
% Forward filtering
for a = 1 : m - 1
    out = input(a, :);
    for idx = 1 : nSections
        y2 = out - states{idx}.z1_A - states{idx}.z2_A;
        misc{idx}.y2(a, :) = y2;
        out = coeff{idx}.b0(a) .* y2 + coeff{idx}.d1(a) .* states{idx}.z1_A + coeff{idx}.d2(a) .* states{idx}.z2_A;
        states{idx}.z2_A = states{idx}.z2_A + coeff{idx}.c2(a) .* states{idx}.z1_A;
        states{idx}.z1_A = states{idx}.z1_A + coeff{idx}.c1(a) .* y2;
        misc{idx}.z2_A(a + 1, :) = states{idx}.z2_A;
        misc{idx}.z1_A(a + 1, :) = states{idx}.z1_A;
    end
    if a == 1
        tmp2 = out;
    else
        tmp2(a, :) = out;
    end
end
out = input(m, :);
for idx = 1 : nSections
    y2 = out - states{idx}.z1_A - states{idx}.z2_A;
    misc{idx}.y2(m, :) = y2;
    out = coeff{idx}.b0(m) .* y2 + coeff{idx}.d1(m) .* states{idx}.z1_A + coeff{idx}.d2(m) .* states{idx}.z2_A;
end
tmp2(m, :) = out;
if reverseFilter
    tmp2 = flipud(tmp2);
end
sub = tmp2 * dp - target;
mse1 = sum(sub.^2);
%% Derivative
sub_grad = sub * 2;
tmp2_grad = sub_grad * dp';
if reverseFilter
    tmp2_grad = flipud(tmp2_grad);
end
grad = cell(nSections, 1);
for idx = 1 : nSections
    grad{idx}.d2_grad = zeros(m, 1);
    grad{idx}.z2_A_grad = zeros(1, nSigs);
    grad{idx}.d1_grad = zeros(m, 1);
    grad{idx}.z1_A_grad = zeros(1, nSigs);
    grad{idx}.b0_grad = zeros(m, 1);
    grad{idx}.y2_grad = zeros(1, nSigs);
    grad{idx}.c1_grad = zeros(m, 1);
    grad{idx}.c2_grad = zeros(m, 1);
end
input_grad = zeros(m, nSigs);
for a = m : -1 : 2
    gd = tmp2_grad(a, :);
    for idx = nSections : -1 : 1
        grad{idx}.d2_grad(a) = grad{idx}.d2_grad(a) + sum(gd .* misc{idx}.z2_A(a, :));
        grad{idx}.d1_grad(a) = grad{idx}.d1_grad(a) + sum(gd .* misc{idx}.z1_A(a, :));
        grad{idx}.b0_grad(a) = grad{idx}.b0_grad(a) + sum(gd .* misc{idx}.y2(a, :));
        grad{idx}.y2_grad = grad{idx}.y2_grad + gd .* coeff{idx}.b0(a);
        grad{idx}.z2_A_grad = grad{idx}.z2_A_grad + gd .* coeff{idx}.d2(a) - grad{idx}.y2_grad;
        grad{idx}.z1_A_grad = grad{idx}.z1_A_grad + gd .* coeff{idx}.d1(a) - grad{idx}.y2_grad;
        gd = grad{idx}.y2_grad;
        grad{idx}.c1_grad(a - 1) = grad{idx}.c1_grad(a - 1) + sum(grad{idx}.z1_A_grad .* misc{idx}.y2(a - 1, :));
        grad{idx}.y2_grad = grad{idx}.z1_A_grad .* coeff{idx}.c1(a - 1);
        grad{idx}.c2_grad(a - 1) = grad{idx}.c2_grad(a - 1) + sum(grad{idx}.z2_A_grad .* misc{idx}.z1_A(a - 1, :));
        grad{idx}.z1_A_grad = grad{idx}.z1_A_grad + grad{idx}.z2_A_grad .* coeff{idx}.c2(a - 1);
    end
    input_grad(a, :) = input_grad(a, :) + gd;
end
gd = tmp2_grad(1, :);
for idx = nSections : -1 : 1
    grad{idx}.d2_grad(1) = grad{idx}.d2_grad(1) + sum(gd .* misc{idx}.z2_A(1, :));
    grad{idx}.d1_grad(1) = grad{idx}.d1_grad(1) + sum(gd .* misc{idx}.z1_A(1, :));
    grad{idx}.b0_grad(1) = grad{idx}.b0_grad(1) + sum(gd .* misc{idx}.y2(1, :));
    grad{idx}.y2_grad = grad{idx}.y2_grad + gd * coeff{idx}.b0(1);
    grad{idx}.z2_A_grad = grad{idx}.z2_A_grad + gd .* coeff{idx}.d2(1) - grad{idx}.y2_grad;
    grad{idx}.z1_A_grad = grad{idx}.z1_A_grad + gd .* coeff{idx}.d1(1) - grad{idx}.y2_grad;
    gd = grad{idx}.y2_grad;
end
input_grad(1, :) = input_grad(1, :) + gd;
if reverseFilter
    for idx = 1 : nSections
        grad{idx}.b0_grad = flipud(grad{idx}.b0_grad);
        grad{idx}.d1_grad = flipud(grad{idx}.d1_grad);
        grad{idx}.d2_grad = flipud(grad{idx}.d2_grad);
        grad{idx}.c1_grad = flipud(grad{idx}.c1_grad);
        grad{idx}.c2_grad = flipud(grad{idx}.c2_grad);
    end
    input_grad = flipud(input_grad);
end
grad3 = [];
for idx = 1 : nSections
    grad3 = [grad3; grad{idx}.b0_grad; grad{idx}.d1_grad; grad{idx}.d2_grad; grad{idx}.c1_grad; grad{idx}.c2_grad];
end
grad3 = [grad3; input_grad(:)];
for idx = 1 : nSections
    grad3 = [grad3; grad{idx}.z2_A_grad(:); grad{idx}.z1_A_grad(:)];
end
plot(grad3);axis tight