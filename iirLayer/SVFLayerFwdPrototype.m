rng(1)
clear
m = 5;
input = randn(m, 2);
dp = randn(2, 1);
polesOrig = randn(m * 5, 1);
%% My implementation
% disp(gradf(1))
poles = ADNode(polesOrig);
clear tmp2 secondOrderOut
b0 = poles(1 : m);
d1 = poles(m + 1 : m * 2);
d2 = poles(m * 2 + 1 : m * 3);
c1 = poles(m * 3 + 1 : m * 4);
c2 = poles(m * 4 + 1 : m * 5);
% Forward filtering
tmp2 = b0(1) .* input(1, :);
z2_A = 0;
z1_A = c1(1) .* input(1, :);
for a = 2 : m - 1
    x = input(a, :);
    y = x - z1_A - z2_A;
    tmp2(a, :) = b0(a) .* y + d1(a) .* z1_A + d2(a) .* z2_A;
    z2_A = z2_A + c2(a) .* z1_A;
    z1_A = z1_A + c1(a) .* y;
end
x = input(m, :);
y = x - z1_A - z2_A;
tmp2(m, :) = b0(m) .* y + d1(m) .* z1_A + d2(m) .* z2_A;
sub = tmp2 * dp;
mse1 = sum(sub.^2);
grad1 = mse1.backprop(1);
%% My implementation
% disp(gradf(1))
rng(1)
poles = [polesOrig; input(:); 0; 0];
poles = ADNode2('base', 'base', 'base', poles);
clear tmp2 secondOrderOut
b0 = poles(1 : m); b0.changeName('b0')
d1 = poles(m + 1 : m * 2); d1.changeName('d1');
d2 = poles(m * 2 + 1 : m * 3); d2.changeName('d2');
c1 = poles(m * 3 + 1 : m * 4); c1.changeName('c1');
c2 = poles(m * 4 + 1 : m * 5); c2.changeName('c2');
flat = poles(m * 5 + 1 : m * 7);
input = reshape(flat, size(input)); input.changeName('input');
% Forward filtering
z2_A = poles(m * 7 + 1); z2_A.changeName('z2_A');
z1_A = poles(m * 7 + 2); z1_A.changeName('z1_A');
for a = 1 : m - 1
    idx = '(' + string(a) + ')';
    y1 = minus(input(a, :), z1_A, 'y1' + idx, 'input' + idx, 'z1_A' + idx);
    y2 = minus(y1, z2_A, 'y2' + idx, 'y1' + idx, 'z2_A' + idx);
    cof1 = b0(a); cof1.changeName('b0' + idx);
    firTerm1 = times(cof1, y2, 'firTerm1' + idx, 'b0' + idx, 'y2' + idx);
    cof2 = d1(a); cof2.changeName('d1' + idx);
    firTerm2 = times(cof2, z1_A, 'firTerm2' + idx, 'd1' + idx, 'z1_A' + idx);
    cof3 = d2(a); cof3.changeName('d2' + idx);
    firTerm3 = times(cof3, z2_A, 'firTerm3' + idx, 'd2' + idx, 'z2_A' + idx);
    sum1_2FirTerm = plus(firTerm1, firTerm2, 'sum1_2FirTerm' + idx, 'firTerm1' + idx, 'firTerm2' + idx);
    out = plus(firTerm3, sum1_2FirTerm, 'out' + idx, 'firTerm3' + idx, 'sum1_2FirTerm' + idx);
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
y1 = minus(input(m, :), z1_A, 'y1' + idx, 'input' + idx, 'z1_A' + idx);
y2 = minus(y1, z2_A, 'y2' + idx, 'y1' + idx, 'z2_A' + idx);
cof1 = b0(m); cof1.changeName('b0' + idx);
firTerm1 = times(cof1, y2, 'firTerm1' + idx, 'b0' + idx, 'y2' + idx);
cof2 = d1(m); cof2.changeName('d1' + idx);
firTerm2 = times(cof2, z1_A, 'firTerm2' + idx, 'd1' + idx, 'z1_A' + idx);
cof3 = d2(m); cof3.changeName('d2' + idx);
firTerm3 = times(cof3, z2_A, 'firTerm3' + idx, 'd2' + idx, 'z2_A' + idx);
sum1_2FirTerm = plus(firTerm1, firTerm2, 'sum1_2FirTerm' + idx, 'firTerm1' + idx, 'firTerm2' + idx);
out = plus(firTerm3, sum1_2FirTerm, 'out' + idx, 'firTerm3' + idx, 'sum1_2FirTerm' + idx);
tmp2(m, :) = out;
sub = tmp2 * dp;
mse2 = sum(sub.^2);
grad2 = mse2.backprop();
% disp(grad1(1))