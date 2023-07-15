rng(1)
input = randn(7, 3);
target = randn(7, 3);

weights = randn(17, 17);

mtx = AutoDiff(weights);
mtx1 = mtx(1 : 7, 1 : 7);
mtx2 = mtx(8 : 10, 8 : 10);
mtx3 = mtx(11 : end, 11 : end);

out = mtx1 * (input * mtx2) + mtx3 * input;

sse = abs(out - target) .^ 2;
flatSse = sum(sse(:));

fval1 = getvalue(flatSse);
grad1 = full(getderivs(flatSse));
grad1 = reshape(grad1, [17, 17]);
%% Matlab sym
% mtxS = sym('mtxS', [17 * 17, 1]);
% 
% mtx = reshape(mtxS, [17, 17]);
% 
% mtx1 = mtx(1 : 7, 1 : 7);
% mtx2 = mtx(8 : 10, 8 : 10);
% mtx3 = mtx(11 : end, 11 : end);
% 
% out = mtx1 * (input * mtx2) + mtx3 * input;
% 
% sse = abs(out - target) .^ 2;
% flatSse = sum(sse(:));
% 
% fh = matlabFunction(flatSse, jacobian(flatSse, mtxS), 'vars', {mtxS});
% 
% [fval2, grad2] = fh(weights(:));
% grad2 = reshape(grad2, [17, 17]);
%% ADNode
mtx = ADNode(weights);

mtx1 = mtx(1 : 7, 1 : 7);
mtx2 = mtx(8 : 10, 8 : 10);
mtx3 = mtx(11 : end, 11 : end);

out = mtx1 * (input * mtx2) + mtx3 * input;

sse = abs(out - target) .^ 2;
flatSse = sum(sse(:));

fval3 = flatSse.value;
grad3 = flatSse.backprop(1);