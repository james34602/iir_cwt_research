function [h_symb, A] = svf_tfMNA(order)
[A, mtxDim] = generate_matrix(order);
Z = inv(A);
h_symb = simplify(Z(mtxDim - 1, mtxDim));
end
function [A, mtxDim] = generate_matrix(order)
N = 3 * order;  % Matrix will be (N+3)x(N+3)
mtxDim = N + 3;
A = sym(zeros(mtxDim, mtxDim));
for i = 1:order
    name = ['d', num2str(i-1)];
    eval(['syms ', name]);  % -di
    name = ['c', num2str(i)];
    eval(['syms ', name]);  % -ci
end
name = ['d', num2str(order)];
eval(['syms ', name]);  % -di
syms z
A(eye(size(A))==1) = 1;
X = 3*(1:order);
Y = 3*(1:order)-1;
for i = 1 : order
    A(X(i), Y(i)) = -1;
    A(X(i) - 2, end) = -1;
    A(X(i), N + 3 - 1) = 1;
    A(X(i) - 1, Y(i) + 1) = -1/z;
    A(X(i), Y(i) + 2) = eval(['-d',num2str(i)]);
    if i < order
        A(X(i), Y(i) + 3) = eval(['-c',num2str(i + 1)]);
    end
end
A(X(order) + 1, end) = -1;
A(N + 3 - 1, 1) = -d0;
A(N + 3 - 1, 2) = -c1;
end