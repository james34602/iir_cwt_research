rng(1)
clear
m = 50;
nSigs = 9;
masterInput = randn(m, nSigs);
dp = randn(nSigs, 1);
polesOrig = rand(m * 5, 1);
reverseFilter = 1;
%% My implementation
% disp(gradf(1))
rng(1)
initialValueZ2 = randn(1, nSigs) * 0.1;
initialValueZ1 = randn(1, nSigs) * 0.1;
poles = [polesOrig; masterInput(:); initialValueZ2(:); initialValueZ1(:)];
poles = ADNode2('base', 'base', 'base', poles);
clear tmp2 secondOrderOut
b0 = poles(1 : m); b0.changeName('b0')
d1 = poles(m + 1 : m * 2); d1.changeName('d1');
d2 = poles(m * 2 + 1 : m * 3); d2.changeName('d2');
c1 = poles(m * 3 + 1 : m * 4); c1.changeName('c1');
c2 = poles(m * 4 + 1 : m * 5); c2.changeName('c2');
flat = poles(m * 5 + 1 : m * 5 + m * nSigs);
input = reshape(flat, size(masterInput)); input.changeName('input');
% Forward filtering
z2_A = poles(m * 5 + m * nSigs + 1 : m * 5 + m * nSigs + nSigs);z2_A = reshape(z2_A, fliplr(size(z2_A.value)));z2_A.changeName('z2_A(0)');
z1_A = poles(m * 5 + m * nSigs + nSigs + 1 : m * 5 + m * nSigs + 2 * nSigs);z1_A = reshape(z1_A, fliplr(size(z1_A.value)));z1_A.changeName('z1_A(0)');
if reverseFilter
    for a = m : -1 : 2
        idx = '(' + string(a) + ')';
        cof0 = input(a, :); cof0.changeName('input' + idx);
        y1 = minus(cof0, z1_A, 'y1' + idx, 'input' + idx, 'z1_A' + idx);
        y2 = minus(y1, z2_A, 'y2' + idx, 'y1' + idx, 'z2_A' + idx);
        cof1 = b0(a); cof1.changeName('b0' + idx);
        firTerm1 = times(cof1, y2, 'firTerm1' + idx, 'b0' + idx, 'y2' + idx);
        cof2 = d1(a); cof2.changeName('d1' + idx);
        firTerm2 = times(cof2, z1_A, 'firTerm2' + idx, 'd1' + idx, 'z1_A' + idx);
        cof3 = d2(a); cof3.changeName('d2' + idx);
        firTerm3 = times(cof3, z2_A, 'firTerm3' + idx, 'd2' + idx, 'z2_A' + idx);
        sum1_2FirTerm = plus(firTerm1, firTerm2, 'sum1_2FirTerm' + idx, 'firTerm1' + idx, 'firTerm2' + idx);
        out = plus(firTerm3, sum1_2FirTerm, 'tmp2(' + string(m - a + 1) + ', :)', 'firTerm3' + idx, 'sum1_2FirTerm' + idx);
        if a == m
            tmp2 = out;
        else
            tmp2(m - a + 1, :) = out;
        end
        cof4 = c2(a); cof4.changeName('c2' + idx);
        wd1 = times(cof4, z1_A, 'wd1' + idx, 'c2' + idx, 'z1_A' + idx);
        z2_A = plus(z2_A, wd1, 'z2_A' + idx, 'z2_A' + idx, 'wd1' + idx);
        cof5 = c1(a); cof5.changeName('c1' + idx);
        wd2 = times(cof5, y2, 'wd2' + idx, 'c1' + idx, 'y2' + idx);
        z1_A = plus(z1_A, wd2, 'z1_A' + idx, 'z1_A' + idx, 'wd2' + idx);
    end
    idx = '(' + string(1) + ')';
    cof0 = input(1, :); cof0.changeName('input' + idx);
    y1 = minus(cof0, z1_A, 'y1' + idx, 'input' + idx, 'z1_A' + idx);
    y2 = minus(y1, z2_A, 'y2' + idx, 'y1' + idx, 'z2_A' + idx);
    cof1 = b0(1); cof1.changeName('b0' + idx);
    firTerm1 = times(cof1, y2, 'firTerm1' + idx, 'b0' + idx, 'y2' + idx);
    cof2 = d1(1); cof2.changeName('d1' + idx);
    firTerm2 = times(cof2, z1_A, 'firTerm2' + idx, 'd1' + idx, 'z1_A' + idx);
    cof3 = d2(1); cof3.changeName('d2' + idx);
    firTerm3 = times(cof3, z2_A, 'firTerm3' + idx, 'd2' + idx, 'z2_A' + idx);
    sum1_2FirTerm = plus(firTerm1, firTerm2, 'sum1_2FirTerm' + idx, 'firTerm1' + idx, 'firTerm2' + idx);
    out = plus(firTerm3, sum1_2FirTerm, 'tmp2(' + string(m) + ', :)', 'firTerm3' + idx, 'sum1_2FirTerm' + idx);
else
    for a = 1 : m - 1
        idx = '(' + string(a) + ')';
        cof0 = input(a, :); cof0.changeName('input' + idx);
        y1 = minus(cof0, z1_A, 'y1' + idx, 'input' + idx, 'z1_A' + idx);
        y2 = minus(y1, z2_A, 'y2' + idx, 'y1' + idx, 'z2_A' + idx);
        cof1 = b0(a); cof1.changeName('b0' + idx);
        firTerm1 = times(cof1, y2, 'firTerm1' + idx, 'b0' + idx, 'y2' + idx);
        cof2 = d1(a); cof2.changeName('d1' + idx);
        firTerm2 = times(cof2, z1_A, 'firTerm2' + idx, 'd1' + idx, 'z1_A' + idx);
        cof3 = d2(a); cof3.changeName('d2' + idx);
        firTerm3 = times(cof3, z2_A, 'firTerm3' + idx, 'd2' + idx, 'z2_A' + idx);
        sum1_2FirTerm = plus(firTerm1, firTerm2, 'sum1_2FirTerm' + idx, 'firTerm1' + idx, 'firTerm2' + idx);
        out = plus(firTerm3, sum1_2FirTerm, 'tmp2(' + string(a) + ', :)', 'firTerm3' + idx, 'sum1_2FirTerm' + idx);
        if a == 1
            tmp2 = out;
        else
            tmp2(a, :) = out;
        end
        cof4 = c2(a); cof4.changeName('c2' + idx);
        wd1 = times(cof4, z1_A, 'wd1' + idx, 'c2' + idx, 'z1_A' + idx);
        z2_A = plus(z2_A, wd1, 'z2_A' + idx, 'z2_A' + idx, 'wd1' + idx);
        cof5 = c1(a); cof5.changeName('c1' + idx);
        wd2 = times(cof5, y2, 'wd2' + idx, 'c1' + idx, 'y2' + idx);
        z1_A = plus(z1_A, wd2, 'z1_A' + idx, 'z1_A' + idx, 'wd2' + idx);
    end
    idx = '(' + string(m) + ')';
    cof0 = input(m, :); cof0.changeName('input' + idx);
    y1 = minus(cof0, z1_A, 'y1' + idx, 'input' + idx, 'z1_A' + idx);
    y2 = minus(y1, z2_A, 'y2' + idx, 'y1' + idx, 'z2_A' + idx);
    cof1 = b0(m); cof1.changeName('b0' + idx);
    firTerm1 = times(cof1, y2, 'firTerm1' + idx, 'b0' + idx, 'y2' + idx);
    cof2 = d1(m); cof2.changeName('d1' + idx);
    firTerm2 = times(cof2, z1_A, 'firTerm2' + idx, 'd1' + idx, 'z1_A' + idx);
    cof3 = d2(m); cof3.changeName('d2' + idx);
    firTerm3 = times(cof3, z2_A, 'firTerm3' + idx, 'd2' + idx, 'z2_A' + idx);
    sum1_2FirTerm = plus(firTerm1, firTerm2, 'sum1_2FirTerm' + idx, 'firTerm1' + idx, 'firTerm2' + idx);
    out = plus(firTerm3, sum1_2FirTerm, 'tmp2(' + string(m) + ', :)', 'firTerm3' + idx, 'sum1_2FirTerm' + idx);
end
tmp2(m, :) = out;
sub = flipud(tmp2) * dp;
mse2 = sum(sub.^2);
grad2 = mse2.backprop();
% grad2 = grad2(1 : m * 5 + m * nSigs);
% disp(grad1(1))
%%
clear y1 y2 firTerm1 firTerm2 firTerm3 sum1_2FirTerm out tmp2 wd1 wd2
poles = [polesOrig; masterInput(:)];
b0 = poles(1 : m);
d1 = poles(m + 1 : m * 2);
d2 = poles(m * 2 + 1 : m * 3);
c1 = poles(m * 3 + 1 : m * 4);
c2 = poles(m * 4 + 1 : m * 5);
% Forward filtering
if reverseFilter
    masterInput = flipud(masterInput);
    b0 = flipud(b0);
    d1 = flipud(d1);
    d2 = flipud(d2);
    c1 = flipud(c1);
    c2 = flipud(c2);
end
z2_A = initialValueZ2; z1_A = initialValueZ1;
tmp2 = zeros(m, nSigs);
for a = 1 : m - 1
    y2(a, :) = masterInput(a, :) - z1_A(a, :) - z2_A(a, :);
    tmp2(a, :) = b0(a) .* y2(a, :) + d1(a) .* z1_A(a, :) + d2(a) .* z2_A(a, :);
    z2_A(a + 1, :) = z2_A(a, :) + c2(a) .* z1_A(a, :);
    z1_A(a + 1, :) = z1_A(a, :) + c1(a) .* y2(a, :);
end
y2(m, :) = masterInput(m, :) - z1_A(m, :) - z2_A(m, :);
tmp2(m, :) = b0(m) .* y2(m, :) + d1(m) .* z1_A(m, :) + d2(m) .* z2_A(m, :);
sub = tmp2 * dp;
mse1 = sum(sub.^2);
%% Grad
sub_grad = sub * 2;
tmp2_grad = sub_grad * dp';
d2_grad = zeros(m, 1);
z2_A_grad = zeros(1, nSigs);
d1_grad = zeros(m, 1);
z1_A_grad = zeros(1, nSigs);
b0_grad = zeros(m, 1);
y2_grad = zeros(1, nSigs);
c1_grad = zeros(m, 1);
c2_grad = zeros(m, 1);
input_grad = zeros(m, nSigs);
for a = m : -1 : 2
    d2_grad(a) = d2_grad(a) + sum(tmp2_grad(a, :) .* z2_A(a, :));
    d1_grad(a) = d1_grad(a) + sum(tmp2_grad(a, :) .* z1_A(a, :));
    b0_grad(a) = b0_grad(a) + sum(tmp2_grad(a, :) .* y2(a, :));
    y2_grad = y2_grad + tmp2_grad(a, :) .* b0(a);
    z2_A_grad = z2_A_grad + tmp2_grad(a, :) .* d2(a) - y2_grad;
    z1_A_grad = z1_A_grad + tmp2_grad(a, :) .* d1(a) - y2_grad;
    input_grad(a, :) = input_grad(a, :) + y2_grad;
    c1_grad(a - 1) = c1_grad(a - 1) + sum(z1_A_grad .* y2(a - 1, :));
    y2_grad = z1_A_grad .* c1(a - 1);
    c2_grad(a - 1) = c2_grad(a - 1) + sum(z2_A_grad .* z1_A(a - 1, :));
    z1_A_grad = z1_A_grad + z2_A_grad .* c2(a - 1);
end
d2_grad(1) = d2_grad(1) + sum(tmp2_grad(1, :) .* initialValueZ2);
d1_grad(1) = d1_grad(1) + sum(tmp2_grad(1, :) .* initialValueZ1);
b0_grad(1) = b0_grad(1) + sum(tmp2_grad(1, :) .* y2(1, :));
y2_grad = y2_grad + tmp2_grad(1, :) * b0(1);
z2_A_grad = z2_A_grad + tmp2_grad(1, :) .* d2(1) - y2_grad;
z1_A_grad = z1_A_grad + tmp2_grad(1, :) .* d1(1) - y2_grad;
input_grad(1, :) = input_grad(1, :) + y2_grad;
if reverseFilter
    b0_grad = flipud(b0_grad);
    d1_grad = flipud(d1_grad);
    d2_grad = flipud(d2_grad);
    c1_grad = flipud(c1_grad);
    c2_grad = flipud(c2_grad);
    input_grad = flipud(input_grad);
end
grad3 = [b0_grad; d1_grad; d2_grad; c1_grad; c2_grad; input_grad(:); z2_A_grad'; z1_A_grad'];
plot(grad2);hold on;plot(grad3);hold off;axis tight
% disp(sum(abs(grad2 - grad3)))