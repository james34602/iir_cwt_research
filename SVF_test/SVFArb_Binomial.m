halfLen = 4097;
halfWndLen = halfLen - 1;
digw = linspace(0, pi - pi / halfWndLen, halfWndLen);
digw(halfWndLen) = pi - pi / halfWndLen;
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
syms c1 c2 c3 c4 c5 c6 c7 c8 d0 d1 d2 d3 d4 d5 d6 d7 d8
% syms b0 b1 b2 b3 b4 b5 b6 b7 b8 a1 a2 a3 a4 a5 a6 a7 a8
z1 = cplxFreq.^-1;
z2 = cplxFreq.^-2;
z3 = cplxFreq.^-3;
z4 = cplxFreq.^-4;
z5 = cplxFreq.^-5;
z6 = cplxFreq.^-6;
z7 = cplxFreq.^-7;
z8 = cplxFreq.^-8;
%% Fourth order
order = 20;
mtx = tf2svfMtx(order);
[b,a] = ellip(4,10,50,0.1,'low');
% [b, a] = cheby1(2, 10, [0.1, 0.6], 'bandpass');
b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);b4 = b(5);
a1 = a(2);a2 = a(3);a3 = a(4);a4 = a(5);
c1 = a1 + 4;
c2 = (6 + 3 * a1 + a2) / c1;
c3 = (4 + 3 * a1 + 2 * a2 + a3) / (c1 * c2);
c4 = (1 + a1 + a2 + a3 + a4) / (c1 * c2 * c3);
d0 = b0;
d1 = (4 * b0 + b1) / c1;
d2 = (6 * b0 + 3 * b1 + b2) / (c1 * c2);
d3 = (4 * b0 + 3 * b1 + 2 * b2 + b3) / (c1 * c2 * c3);
d4 = (b0 + b1 + b2 + b3 + b4) / (c1 * c2 * c3 * c4);
unnormalizedc = mtx * a';
for c = 1 : order
    unnormalizedc(c + 1) = unnormalizedc(c + 1) / prod(unnormalizedc(1 : c));
end
unnormalizedd = mtx * b';
for c = 1 : order + 1
    unnormalizedd(c) = unnormalizedd(c) / prod(unnormalizedc(1 : c));
end
%% Eighth order
order = 8;
[mtx, mtxInv] = tf2svfMtx(order);
[b, a] = ellip(4, 10, 50, [0.3, 0.6], 'bandpass');
unnormalizedc = mtx * a';
for c = 1 : order
    unnormalizedc(c + 1) = unnormalizedc(c + 1) / prod(unnormalizedc(1 : c));
end
unnormalizedd = mtx * b';
for c = 1 : order + 1
    unnormalizedd(c) = unnormalizedd(c) / prod(unnormalizedc(1 : c));
end
for c = 1 : order
    eval("c" + string(c) + "=unnormalizedc("+string(c+1)+")");
    eval("d" + string(c - 1) + "=unnormalizedd("+string(c)+")");
end
eval("d" + string(order) + "=unnormalizedd("+string(order+1)+")");
tfEig = (d0 - 8*d0*z1 + 28*d0*z2 - 56*d0*z3 + 70*d0*z4 - 56*d0*z5 + 28*d0*z6 - 8*d0*z7 + d0*z8 + c1*d1*z1 - 7*c1*d1*z2 + 21*c1*d1*z3 - 35*c1*d1*z4 + 35*c1*d1*z5 - 21*c1*d1*z6 + 7*c1*d1*z7 - c1*d1*z8 + c1*c2*d2*z2 - 6*c1*c2*d2*z3 + 15*c1*c2*d2*z4 - 20*c1*c2*d2*z5 + 15*c1*c2*d2*z6 - 6*c1*c2*d2*z7 + c1*c2*d2*z8 + c1*c2*c3*d3*z3 - 5*c1*c2*c3*d3*z4 + 10*c1*c2*c3*d3*z5 - 10*c1*c2*c3*d3*z6 + 5*c1*c2*c3*d3*z7 - c1*c2*c3*d3*z8 + c1*c2*c3*c4*d4*z4 - 4*c1*c2*c3*c4*d4*z5 + 6*c1*c2*c3*c4*d4*z6 - 4*c1*c2*c3*c4*d4*z7 + c1*c2*c3*c4*d4*z8 + c1*c2*c3*c4*c5*d5*z5 - 3*c1*c2*c3*c4*c5*d5*z6 + 3*c1*c2*c3*c4*c5*d5*z7 - c1*c2*c3*c4*c5*d5*z8 + c1*c2*c3*c4*c5*c6*d6*z6 - 2*c1*c2*c3*c4*c5*c6*d6*z7 + c1*c2*c3*c4*c5*c6*d6*z8 + c1*c2*c3*c4*c5*c6*c7*d7*z7 - c1*c2*c3*c4*c5*c6*c7*d7*z8 + c1*c2*c3*c4*c5*c6*c7*c8*d8*z8)./(c1*z1 - 8*z1 - 7*c1*z2 + 21*c1*z3 - 35*c1*z4 + 35*c1*z5 - 21*c1*z6 + 7*c1*z7 - c1*z8 + 28*z2 - 56*z3 + 70*z4 - 56*z5 + 28*z6 - 8*z7 + z8 + c1*c2*z2 - 6*c1*c2*z3 + 15*c1*c2*z4 - 20*c1*c2*z5 + 15*c1*c2*z6 - 6*c1*c2*z7 + c1*c2*z8 + c1*c2*c3*z3 - 5*c1*c2*c3*z4 + 10*c1*c2*c3*z5 - 10*c1*c2*c3*z6 + 5*c1*c2*c3*z7 - c1*c2*c3*z8 + c1*c2*c3*c4*z4 - 4*c1*c2*c3*c4*z5 + 6*c1*c2*c3*c4*z6 - 4*c1*c2*c3*c4*z7 + c1*c2*c3*c4*z8 + c1*c2*c3*c4*c5*z5 - 3*c1*c2*c3*c4*c5*z6 + 3*c1*c2*c3*c4*c5*z7 - c1*c2*c3*c4*c5*z8 + c1*c2*c3*c4*c5*c6*z6 - 2*c1*c2*c3*c4*c5*c6*z7 + c1*c2*c3*c4*c5*c6*z8 + c1*c2*c3*c4*c5*c6*c7*z7 - c1*c2*c3*c4*c5*c6*c7*z8 + c1*c2*c3*c4*c5*c6*c7*c8*z8 + 1);
imagesc(log(mtx))
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