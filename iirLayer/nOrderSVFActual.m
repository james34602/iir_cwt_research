rng(1)
clear
m = 4;
nSigs = 3;
input = randn(m, nSigs);
dp = randn(nSigs, 1);
nSections = 2;
ltv = 0;
if ltv == 1
    coefNum = m;
else
    coefNum = 1;
end
parallelFilter = 0;
reverseFilter = 1;
statesNoGrad = 0;
target = randn(m, 1);
%% My implementation
b0 = rand(coefNum, nSections);
d1 = rand(coefNum, nSections);
d2 = rand(coefNum, nSections);
c1 = rand(coefNum, nSections);
c2 = rand(coefNum, nSections);
if statesNoGrad
    initialValueZ1 = zeros(nSigs, nSections);
    initialValueZ2 = zeros(nSigs, nSections);
else
    initialValueZ1 = rand(nSigs, nSections);
    initialValueZ2 = rand(nSigs, nSections);
end
poles = [input(:); b0(:); d1(:); d2(:); c1(:); c2(:); initialValueZ1(:); initialValueZ2(:)];
input = reshape(poles(1 : m * nSigs), [m, nSigs]);
b0 = reshape(poles(m * nSigs + coefNum * nSections * 0 + 1 : m * nSigs + coefNum * nSections * 1), [coefNum, nSections]);
d1 = reshape(poles(m * nSigs + coefNum * nSections * 1 + 1 : m * nSigs + coefNum * nSections * 2), [coefNum, nSections]);
d2 = reshape(poles(m * nSigs + coefNum * nSections * 2 + 1 : m * nSigs + coefNum * nSections * 3), [coefNum, nSections]);
c1 = reshape(poles(m * nSigs + coefNum * nSections * 3 + 1 : m * nSigs + coefNum * nSections * 4), [coefNum, nSections]);
c2 = reshape(poles(m * nSigs + coefNum * nSections * 4 + 1 : m * nSigs + coefNum * nSections * 5), [coefNum, nSections]);
initialValueZ1 = reshape(poles(m * nSigs + coefNum * nSections * 5 + 1 : m * nSigs + coefNum * nSections * 5 + nSigs * nSections), [nSigs, nSections]);
initialValueZ2 = reshape(poles(m * nSigs + coefNum * nSections * 5 + nSigs * nSections + 1 : m * nSigs + coefNum * nSections * 5 + nSigs * nSections * 2), [nSigs, nSections]);
z2_A = zeros(m, nSigs, nSections);z1_A = zeros(m, nSigs, nSections);
z1_A(1, :, :) = initialValueZ1;
z2_A(1, :, :) = initialValueZ2;
y2 = zeros(m, nSigs, nSections);
if reverseFilter
    input = flipud(input);b0 = flipud(b0);d1 = flipud(d1);d2 = flipud(d2);c1 = flipud(c1);c2 = flipud(c2);
end
% Forward filtering
tmp2 = zeros(m, nSigs);
if parallelFilter == 1
    out = zeros(nSections, nSigs);
    for a = 1 : m - 1
        if ltv == 1
            coefIdx = a;
        else
            coefIdx = 1;
        end
        in = input(a, :);
        for idx = 1 : nSections
            y2(a, :, idx) = in - z1_A(a, :, idx) - z2_A(a, :, idx);
            out(idx, :) = b0(coefIdx, idx) .* y2(a, :, idx) + d1(coefIdx, idx) .* z1_A(a, :, idx) + d2(coefIdx, idx) .* z2_A(a, :, idx);
            z2_A(a + 1, :, idx) = z2_A(a, :, idx) + c2(coefIdx, idx) .* z1_A(a, :, idx);
            z1_A(a + 1, :, idx) = z1_A(a, :, idx) + c1(coefIdx, idx) .* y2(a, :, idx);
        end
        tmp2(a, :) = sum(out, 1);
    end
    if ltv == 1
        coefIdx = m;
    else
        coefIdx = 1;
    end
    in = input(m, :);
    for idx = 1 : nSections
        y2(m, :, idx) = in - z1_A(m, :, idx) - z2_A(m, :, idx);
        out(idx, :) = b0(coefIdx, idx) .* y2(m, :, idx) + d1(coefIdx, idx) .* z1_A(m, :, idx) + d2(coefIdx, idx) .* z2_A(m, :, idx);
    end
    tmp2(m, :) = sum(out, 1);
else
    for a = 1 : m - 1
        if ltv == 1
            coefIdx = a;
        else
            coefIdx = 1;
        end
        in = input(a, :);
        for idx = 1 : nSections
            y2(a, :, idx) = in - z1_A(a, :, idx) - z2_A(a, :, idx);
            in = b0(coefIdx, idx) .* y2(a, :, idx) + d1(coefIdx, idx) .* z1_A(a, :, idx) + d2(coefIdx, idx) .* z2_A(a, :, idx);
            z2_A(a + 1, :, idx) = z2_A(a, :, idx) + c2(coefIdx, idx) .* z1_A(a, :, idx);
            z1_A(a + 1, :, idx) = z1_A(a, :, idx) + c1(coefIdx, idx) .* y2(a, :, idx);
        end
        tmp2(a, :) = in;
    end
    if ltv == 1
        coefIdx = m;
    else
        coefIdx = 1;
    end
    in = input(m, :);
    for idx = 1 : nSections
        y2(m, :, idx) = in - z1_A(m, :, idx) - z2_A(m, :, idx);
        in = b0(coefIdx, idx) .* y2(m, :, idx) + d1(coefIdx, idx) .* z1_A(m, :, idx) + d2(coefIdx, idx) .* z2_A(m, :, idx);
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
b0_grad = zeros(coefNum, nSections);
d1_grad = zeros(coefNum, nSections);
d2_grad = zeros(coefNum, nSections);
c1_grad = zeros(coefNum, nSections);
c2_grad = zeros(coefNum, nSections);
z2_A_grad = zeros(nSections, nSigs);
z1_A_grad = zeros(nSections, nSigs);
y2_grad = zeros(nSections, nSigs);
section_grad = zeros(nSections, nSigs);
input_grad = zeros(m, nSigs);
for a = m : -1 : 2
    if ltv == 1
        coefIdx1 = a;
        coefIdx2 = a - 1;
    else
        coefIdx1 = 1;
        coefIdx2 = 1;
    end
    gd = tmp2_grad(a, :);
    for idx = nSections : -1 : 1
        d2_grad(coefIdx1, idx) = d2_grad(coefIdx1, idx) + sum(gd .* z2_A(a, :, idx));
        d1_grad(coefIdx1, idx) = d1_grad(coefIdx1, idx) + sum(gd .* z1_A(a, :, idx));
        b0_grad(coefIdx1, idx) = b0_grad(coefIdx1, idx) + sum(gd .* y2(a, :, idx));
        section_grad(idx, :) = y2_grad(idx, :) + gd .* b0(coefIdx1, idx);
        z2_A_grad(idx, :) = z2_A_grad(idx, :) + gd .* d2(coefIdx1, idx) - section_grad(idx, :);
        z1_A_grad(idx, :) = z1_A_grad(idx, :) + gd .* d1(coefIdx1, idx) - section_grad(idx, :);
        if parallelFilter == 0
            gd = section_grad(idx, :);
        end
        c1_grad(coefIdx2, idx) = c1_grad(coefIdx2, idx) + sum(z1_A_grad(idx, :) .* y2(a - 1, :, idx));
        y2_grad(idx, :) = z1_A_grad(idx, :) .* c1(coefIdx2, idx);
        c2_grad(coefIdx2, idx) = c2_grad(coefIdx2, idx) + sum(z2_A_grad(idx, :) .* z1_A(a - 1, :, idx));
        z1_A_grad(idx, :) = z1_A_grad(idx, :) + z2_A_grad(idx, :) .* c2(coefIdx2, idx);
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
z2_A_grad = z2_A_grad.';
z1_A_grad = z1_A_grad.';
%% Actual SVF layer end
grad3 = [input_grad(:); b0_grad(:); d1_grad(:); d2_grad(:); c1_grad(:); c2_grad(:)];
if ~statesNoGrad
    grad3 = [grad3; z1_A_grad(:); z2_A_grad(:)];
end
plot(grad3);axis tight