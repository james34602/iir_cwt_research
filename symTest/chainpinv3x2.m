rng(1)
input = randn(3, 1);
target = randn(2, 1);

weights = randn(3, 2);
%% Matlab sym
mtxS = sym('mtx', [3 * 2, 1]);

mtx = reshape(mtxS, [3, 2]);

out = pinv(mtx) * input;

sse = abs(out - target) .^ 2;
flatSse = sum(sse(:));

fh = matlabFunction(flatSse, jacobian(flatSse, mtxS), 'vars', {mtxS});

[fval2, grad2] = fh(weights(:));
grad2 = reshape(grad2, [3, 2]);
%% ADNode
mtx = ADNode(weights);

out = pinv(mtx) * input;

sse = abs(out - target) .^ 2;
flatSse = sum(sse(:));

fval3 = flatSse.value;
grad3 = flatSse.backprop(1);