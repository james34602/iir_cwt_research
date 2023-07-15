rng(1)
clear
m = 4;
x_fft = randn(m, 2);
% x_fft = sym('x_fft', [m, 1]);
poles = sym('poles', [m * 5, 1]);
b0 = poles(1 : m);
d1 = poles(m + 1 : m * 2);
d2 = poles(m * 2 + 1 : m * 3);
c1 = poles(m * 3 + 1 : m * 4);
c2 = poles(m * 4 + 1 : m * 5);
target = randn(m, 1);
dp = randn(2, 1);
% target = sym('target', [m, 1]);

% Forward filtering
tmp2 = b0(1) * x_fft(1, :);
z2_A = 0;
z1_A = c1(1) * x_fft(1, :);
for a = 2 : m
    x = x_fft(a, :);
    y = x - z1_A - z2_A;
    tmp2(a, :) = b0(a) * y + d1(a) * z1_A + d2(a) * z2_A;
    z2_A = z2_A + c2(a) * z1_A;
    z1_A = z1_A + c1(a) * y;
end
% Reverse filtering
secondOrderOut(m, :) = b0(m) * tmp2(m, :);
z2_A = 0;
z1_A = c1(m) * tmp2(m, :);
for a = m - 1 : -1 : 1
    x = tmp2(a, :);
    y = x - z1_A - z2_A;
    secondOrderOut(a, :) = b0(a) * y + d1(a) * z1_A + d2(a) * z2_A;
    z2_A = z2_A + c2(a) * z1_A;
    z1_A = z1_A + c1(a) * y;
end
sub = secondOrderOut * dp - target;
mse = mean(sub.^2);
gradf = jacobian(mse, poles);
fh = matlabFunction(mse,gradf,'vars',{poles});
poles = randn(m * 5, 1);
[mse, gradf] = fh(poles);
%% My implementation
% disp(gradf(1))
poles = ADNode(poles);
clear tmp2 secondOrderOut
b0 = poles(1 : m);
d1 = poles(m + 1 : m * 2);
d2 = poles(m * 2 + 1 : m * 3);
c1 = poles(m * 3 + 1 : m * 4);
c2 = poles(m * 4 + 1 : m * 5);
% Forward filtering
tmp2 = b0(1) .* x_fft(1, :);
z2_A = 0;
z1_A = c1(1) .* x_fft(1, :);
for a = 2 : m - 1
    x = x_fft(a, :);
    y = x - z1_A - z2_A;
    tmp2(a, :) = b0(a) .* y + d1(a) .* z1_A + d2(a) .* z2_A;
    z2_A = z2_A + c2(a) .* z1_A;
    z1_A = z1_A + c1(a) .* y;
end
x = x_fft(m, :);
y = x - z1_A - z2_A;
tmp2(m, :) = b0(m) .* y + d1(m) .* z1_A + d2(m) .* z2_A;
% Reverse filtering
secondOrderOut = b0(m) .* tmp2(m, :);
z2_A = 0;
z1_A = c1(m) .* tmp2(m, :);
for a = m - 1 : -1 : 2
    x = tmp2(a, :);
    y = x - z1_A - z2_A;
    secondOrderOut(m - a + 1, :) = b0(a) .* y + d1(a) .* z1_A + d2(a) .* z2_A;
    z2_A = z2_A + c2(a) .* z1_A;
    z1_A = z1_A + c1(a) .* y;
end
x = tmp2(1, :);
y = x - z1_A - z2_A;
secondOrderOut(m - 1 + 1, :) = b0(1) .* y + d1(1) .* z1_A + d2(1) .* z2_A;
sub = flipud(secondOrderOut) * dp - target;
mse2 = mean(sub.^2);
grad1 = mse2.backprop(1)
plot(gradf);hold on;plot(grad1);hold off;axis tight
% disp(grad1(1))