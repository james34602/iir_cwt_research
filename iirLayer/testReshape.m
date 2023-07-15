rng(1)
clear
m = 24;
x_fft = randn(m, 1);
vec1 = rand(2, 1);
mtx1 = randn(4, 8);
mtx2 = rand(4, 3);
target = randn(1, 32);
x = sym('x', [m, 1]);
dm1 = reshape(x, [m/2, 2]);
df1 = cos(dm1 * vec1);
dm2 = reshape(df1, [3, 4]);
df2 = sin(dm2 * mtx1);
df3 = mtx2 * df2;
df4 = reshape(df3, [1, 32]);
dif = df4 - target;
mse = mean(dif(:).^2);
gradf = jacobian(mse, x);
fh = matlabFunction(mse,gradf,'vars',{x});
poles = randn(m, 1);
[mse, gradf] = fh(poles);
%% My implementation
% disp(gradf(1))
x = ADNode(poles);
dm1 = reshape(x, [m/2, 2]);
df1 = cos(dm1 * vec1);
dm2 = reshape(df1, [3, 4]);
df2 = sin(dm2 * mtx1);
df3 = mtx2 * df2;
df4 = reshape(df3, [1, 32]);
dif = df4 - target;
mse2 = mean(dif(:).^2);
grad1 = mse2.backprop(1)
% disp(grad1(1))