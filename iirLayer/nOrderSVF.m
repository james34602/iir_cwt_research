rng(1)
clear
m = 40;
nSigs = 3;
masterInput = randn(m, nSigs);
dp = randn(nSigs, 1);
nSections = 2;
ltv = 1;
if ltv == 1
    polesOrig = rand(m * 5 * nSections, 1);
else
    polesOrig = rand(5 * nSections, 1);
end
parallelFilter = 1;
reverseFilter = 1;
target = randn(m, 1);
%% My implementation
poles = [polesOrig; masterInput(:)];
for idx = 1 : nSections
    initialValueZ2 = randn(1, nSigs) * 0.1;
    initialValueZ1 = randn(1, nSigs) * 0.1;
    poles = [poles; initialValueZ2(:); initialValueZ1(:);];
end
if ltv == 0
    polesBk = poles(5 * nSections + m * nSigs + 1 : end);
    poles2 = [];
    for idx = 1 : 5 * nSections
        poles2 = [poles2; repmat(poles(idx), m, 1)];
    end
    poles = [poles2; masterInput(:); polesBk];
end
flat = poles(m * 5 * (nSections - 1) + m * 5 + 1 : m * 5 * (nSections - 1) + m * 5 + m * nSigs);
input = reshape(flat, size(masterInput));
z2_A = zeros(m, nSigs, nSections);z1_A = zeros(m, nSigs, nSections);
b0 = zeros(m, nSections);d1 = zeros(m, nSections);d2 = zeros(m, nSections);c1 = zeros(m, nSections);c2 = zeros(m, nSections);
for idx = 1 : nSections
    base = m * 5 * (idx - 1);
    b0(:, idx) = poles(base + 1 : base + m);
    d1(:, idx) = poles(base + m + 1 : base + m * 2);
    d2(:, idx) = poles(base + m * 2 + 1 : base + m * 3);
    c1(:, idx) = poles(base + m * 3 + 1 : base + m * 4);
    c2(:, idx) = poles(base + m * 4 + 1 : base + m * 5);
    base = m * 5 * (nSections - 1) + m * 5 + m * nSigs + nSigs * 2 * (idx - 1);
    z2_A(1, :, idx) = poles(base + 1 : base + nSigs).';
    z1_A(1, :, idx) = poles(base + nSigs + 1 : base + nSigs * 2).';
end
y2 = zeros(m, nSigs, nSections);
if reverseFilter
    input = flipud(input);b0 = flipud(b0);d1 = flipud(d1);d2 = flipud(d2);c1 = flipud(c1);c2 = flipud(c2);
end
% Forward filtering
tmp2 = zeros(m, nSigs);
if parallelFilter == 1
    out = zeros(nSections, nSigs);
    for a = 1 : m - 1
        in = input(a, :);
        for idx = 1 : nSections
            y2(a, :, idx) = in - z1_A(a, :, idx) - z2_A(a, :, idx);
            out(idx, :) = b0(a, idx) .* y2(a, :, idx) + d1(a, idx) .* z1_A(a, :, idx) + d2(a, idx) .* z2_A(a, :, idx);
            z2_A(a + 1, :, idx) = z2_A(a, :, idx) + c2(a, idx) .* z1_A(a, :, idx);
            z1_A(a + 1, :, idx) = z1_A(a, :, idx) + c1(a, idx) .* y2(a, :, idx);
        end
        tmp2(a, :) = sum(out, 1);
    end
    in = input(m, :);
    for idx = 1 : nSections
        y2(m, :, idx) = in - z1_A(m, :, idx) - z2_A(m, :, idx);
        out(idx, :) = b0(m, idx) .* y2(m, :, idx) + d1(m, idx) .* z1_A(m, :, idx) + d2(m, idx) .* z2_A(m, :, idx);
    end
    tmp2(m, :) = sum(out, 1);
else
    for a = 1 : m - 1
        in = input(a, :);
        for idx = 1 : nSections
            y2(a, :, idx) = in - z1_A(a, :, idx) - z2_A(a, :, idx);
            in = b0(a, idx) .* y2(a, :, idx) + d1(a, idx) .* z1_A(a, :, idx) + d2(a, idx) .* z2_A(a, :, idx);
            z2_A(a + 1, :, idx) = z2_A(a, :, idx) + c2(a, idx) .* z1_A(a, :, idx);
            z1_A(a + 1, :, idx) = z1_A(a, :, idx) + c1(a, idx) .* y2(a, :, idx);
        end
        tmp2(a, :) = in;
    end
    in = input(m, :);
    for idx = 1 : nSections
        y2(m, :, idx) = in - z1_A(m, :, idx) - z2_A(m, :, idx);
        in = b0(m, idx) .* y2(m, :, idx) + d1(m, idx) .* z1_A(m, :, idx) + d2(m, idx) .* z2_A(m, :, idx);
    end
    tmp2(m, :) = in;
end
if reverseFilter
    tmp2 = flipud(tmp2);
end
sub = tmp2 * dp - target;
mse1 = sum(sub.^2);
%% Derivative
sub_grad = sub * 2;
tmp2_grad = sub_grad * dp';
%% Actual SVF layer start
if reverseFilter
    tmp2_grad = flipud(tmp2_grad);
end
b0_grad = zeros(m, nSections);
d1_grad = zeros(m, nSections);
d2_grad = zeros(m, nSections);
c1_grad = zeros(m, nSections);
c2_grad = zeros(m, nSections);
z2_A_grad = zeros(nSections, nSigs);
z1_A_grad = zeros(nSections, nSigs);
y2_grad = zeros(nSections, nSigs);
section_grad = zeros(nSections, nSigs);
input_grad = zeros(m, nSigs);
for a = m : -1 : 2
    gd = tmp2_grad(a, :);
    for idx = nSections : -1 : 1
        d2_grad(a, idx) = d2_grad(a, idx) + sum(gd .* z2_A(a, :, idx));
        d1_grad(a, idx) = d1_grad(a, idx) + sum(gd .* z1_A(a, :, idx));
        b0_grad(a, idx) = b0_grad(a, idx) + sum(gd .* y2(a, :, idx));
        section_grad(idx, :) = y2_grad(idx, :) + gd .* b0(a, idx);
        z2_A_grad(idx, :) = z2_A_grad(idx, :) + gd .* d2(a, idx) - section_grad(idx, :);
        z1_A_grad(idx, :) = z1_A_grad(idx, :) + gd .* d1(a, idx) - section_grad(idx, :);
        if parallelFilter == 0
            gd = section_grad(idx, :);
        end
        c1_grad(a - 1, idx) = c1_grad(a - 1, idx) + sum(z1_A_grad(idx, :) .* y2(a - 1, :, idx));
        y2_grad(idx, :) = z1_A_grad(idx, :) .* c1(a - 1, idx);
        c2_grad(a - 1, idx) = c2_grad(a - 1, idx) + sum(z2_A_grad(idx, :) .* z1_A(a - 1, :, idx));
        z1_A_grad(idx, :) = z1_A_grad(idx, :) + z2_A_grad(idx, :) .* c2(a - 1, idx);
    end
    if parallelFilter == 1
        input_grad(a, :) = input_grad(a, :) + sum(section_grad, 1);
    else
        input_grad(a, :) = input_grad(a, :) + section_grad(1, :);
    end
end
gd = tmp2_grad(1, :);
for idx = nSections : -1 : 1
    d2_grad(1, idx) = d2_grad(1, idx) + sum(gd .* z2_A(1, :, idx));
    d1_grad(1, idx) = d1_grad(1, idx) + sum(gd .* z1_A(1, :, idx));
    b0_grad(1, idx) = b0_grad(1, idx) + sum(gd .* y2(1, :, idx));
    section_grad(idx, :) = y2_grad(idx, :) + gd * b0(1, idx);
    z2_A_grad(idx, :) = z2_A_grad(idx, :) + gd .* d2(1, idx) - section_grad(idx, :);
    z1_A_grad(idx, :) = z1_A_grad(idx, :) + gd .* d1(1, idx) - section_grad(idx, :);
    if parallelFilter == 0
        gd = section_grad(idx, :);
    end
end
if parallelFilter == 1
    input_grad(1, :) = input_grad(1, :) + sum(section_grad, 1);
else
    input_grad(1, :) = input_grad(1, :) + section_grad(1, :);
end
if reverseFilter
    b0_grad = flipud(b0_grad);d1_grad = flipud(d1_grad);d2_grad = flipud(d2_grad);c1_grad = flipud(c1_grad);c2_grad = flipud(c2_grad);input_grad = flipud(input_grad);
end
if ltv == 0
    b0_grad = sum(b0_grad, 1);d1_grad = sum(d1_grad, 1);d2_grad = sum(d2_grad, 1);c1_grad = sum(c1_grad, 1);c2_grad = sum(c2_grad, 1);
end
%% Actual SVF layer end
grad3 = [];
for idx = 1 : nSections
    grad3 = [grad3; b0_grad(:, idx); d1_grad(:, idx); d2_grad(:, idx); c1_grad(:, idx); c2_grad(:, idx)];
end
grad3 = [grad3; input_grad(:)];
for idx = 1 : nSections
    grad3 = [grad3; z2_A_grad(idx, :).'; z1_A_grad(idx, :).'];
end
plot(grad3);axis tight