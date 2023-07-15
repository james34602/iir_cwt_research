rng(1)
clear
m = 24;
n = 16;
q = 40;
mtx1 = randn(8, 10);
target = randn(10, 10);
x = sym('x', [m + n + q, 1]);
A = x(1 : m);
B = x(m + 1 : m + n);
C = x(m + n + 1 : m + n + q);
dm1 = reshape(A, [m/4, 4]);
dm2 = reshape(B, [n/4, 4]);
dm3 = reshape(C, [10, 4]);
df1 = vertcat(dm1, dm2);
df2 = horzcat(df1, dm3);
df4 = df2 * mtx1;
dif = df4 - target;
mse = mean(dif(:).^2);
% gradf = jacobian(mse, x);
% fh = matlabFunction(mse,gradf,'vars',{x});
poles = randn(m + n + q, 1);
% [mse, gradf] = fh(poles);
%% My implementation
% disp(gradf(1))
x = ADNode(poles);
A = x(1 : m);
B = x(m + 1 : m + n);
C = x(m + n + 1 : m + n + q);
dm1 = reshape(A, [m/4, 4]);
dm2 = reshape(B, [n/4, 4]);
dm3 = reshape(C, [10, 4]);
df1 = vertcat(dm1, dm2);
df2 = horzcat(df1, dm3);
df4 = df2 * mtx1;
dif = df4 - target;
mse2 = mean(dif(:).^2);
grad1 = mse2.backprop(1)
% disp(grad1(1))