function [H, mtx, mtxInv] = svf_tfDerived(order)
[mtx, mtxInv] = tf2svfMtx(order);
 aa=circshift(diag(order:-1:0), 1, 1);
 mtx2 = expm(aa);
 mtxInv2 = expm(-aa);
mtxInvSym = sym(mtxInv);
v=['[' sprintf('c%d;',0:order)];v(end)=']';
c=str2sym(v);
for ord = order : -1 : 1
    c(ord + 1) = c(ord + 1) * prod(c(1 : ord));
end
a_hat = mtxInvSym * c;
a_hat = simplify(a_hat/c(1));
v=['[' sprintf('d%d;',0:order)];v(end)=']';
d=str2sym(v);
c2 = [1; sym('c', [order, 1])];
for ord = order + 1 : -1 : 1
    d(ord) = d(ord) * prod(c2(1 : ord));
end
b_hat = mtxInvSym * d;
b_hat = simplify(b_hat);
syms z
numerator = 0;
denominator = 0;
for k = 0:order
    numerator = numerator + b_hat(k+1) * z^(-k);
    denominator = denominator + a_hat(k+1) * z^(-k);
end
H = numerator / denominator;
H = simplify(H);
end
function [mtx, mtxInv] = tf2svfMtx(order)
mtx = zeros(order + 1, order + 1);
for c = 0 : order - 1
    for a = 1 : order
        if c + a <= order
            mtx(order - c + 1, order - a + 1 - c) = nchoosek(c + a,c);
        end
        mtx(order - c + 1, order + 1 - c) = 1;
    end
end
mtx(1, 1) = 1;
mtxInv = zeros(order + 1, order + 1);
for a = 1 : order + 1
    sgn = (-1) ^ (a - 1);
    for b = 1 : order - a + 2
        mtxInv((a - 1) + b, b) = mtx((a - 1) + b, b) * sgn;
    end
end
end