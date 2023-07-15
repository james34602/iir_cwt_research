rng(1)
input = randn(2, 1);
target = randn(3, 1);

weights = randn(2, 3);
%% Matlab sym
mtxS = sym('mtx', [2 * 3, 1]);

mtx = reshape(mtxS, [2, 3]);

out = pinv(mtx) * input;

sse = abs(out - target) .^ 2;
flatSse = sum(sse(:));

fh = matlabFunction(flatSse, jacobian(flatSse, mtxS), 'vars', {mtxS});

[fval2, grad2] = fh(weights(:));
grad2 = reshape(grad2, [2, 3]);
%% ADNode
mtx = ADNode(weights);

out = pinv(mtx) * input;

sse = abs(out - target) .^ 2;
flatSse = sum(sse(:));

fval3 = flatSse.value;
grad3 = flatSse.backprop(1);