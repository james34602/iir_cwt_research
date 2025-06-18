halfLen = 4097;
halfWndLen = halfLen - 1;
digw = linspace(0, pi - pi / halfWndLen, halfWndLen);
digw(halfWndLen) = pi - pi / halfWndLen;
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
syms z
z1 = z^-1;
z2 = z^-2;
z3 = z^-3;
z4 = z^-4;
z5 = z^-5;
z6 = z^-6;
z7 = z^-7;
z8 = z^-8;
syms c1 c2 c3 c4 c5 c6 c7 c8 d0 d1 d2 d3 d4 d5 d6 d7 d8
term = c1*z1 - 2*z1 - c1*z2 + z2 + c1*c2*z2 + 1;
branch1 = (z2 - 2*z1 + 1) ./ term;
branch2 = (c1 * (z1 - z2)) ./ term;
branch3 = (c1*c2 * z2) ./ term;
tfS = d0 * branch1 + d1 * branch2 + d2 * branch3;
tfT = (d0 - 3*d0*z1 + 3*d0*z2 - d0*z3 + c1*d1*z1 - 2*c1*d1*z2 + c1*d1*z3 + c1*c2*d2*z2 - c1*c2*d2*z3 + c1*c2*c3*d3*z3)./(c1*z1 - 3*z1 - 2*c1*z2 + c1*z3 + 3*z2 - z3 + c1*c2*z2 - c1*c2*z3 + c1*c2*c3*z3 + 1);
tfF = (d0 - 4*d0*z1 + 6*d0*z2 - 4*d0*z3 + d0*z4 + c1*d1*z1 - 3*c1*d1*z2 + 3*c1*d1*z3 - c1*d1*z4 + c1*c2*d2*z2 - 2*c1*c2*d2*z3 + c1*c2*d2*z4 + c1*c2*c3*d3*z3 - c1*c2*c3*d3*z4 + c1*c2*c3*c4*d4*z4)./(c1*z1 - 4*z1 - 3*c1*z2 + 3*c1*z3 - c1*z4 + 6*z2 - 4*z3 + z4 + c1*c2*z2 - 2*c1*c2*z3 + c1*c2*z4 + c1*c2*c3*z3 - c1*c2*c3*z4 + c1*c2*c3*c4*z4 + 1);
tfFive = (c1*d1 - d0 + 5*d0*z - 10*d0*z^2 + 10*d0*z^3 - 5*d0*z^4 + d0*z^5 - c1*c2*d2 - 4*c1*d1*z + 6*c1*d1*z^2 - 4*c1*d1*z^3 + c1*d1*z^4 + 3*c1*c2*d2*z - 3*c1*c2*d2*z^2 + c1*c2*d2*z^3 + c1*c2*c3*d3 - c1*c2*c3*c4*d4 - 2*c1*c2*c3*d3*z + c1*c2*c3*d3*z^2 + c1*c2*c3*c4*c5*d5 + c1*c2*c3*c4*d4*z)/(c1 + 5*z - c1*c2 - 4*c1*z + 6*c1*z^2 - 4*c1*z^3 + c1*z^4 - 10*z^2 + 10*z^3 - 5*z^4 + z^5 + c1*c2*c3 + 3*c1*c2*z - 3*c1*c2*z^2 + c1*c2*z^3 - 2*c1*c2*c3*z + c1*c2*c3*z^2 - c1*c2*c3*c4 + c1*c2*c3*c4*c5 + c1*c2*c3*c4*z - 1);
tfSix = (d0 - c1*d1 - 6*d0*z + 15*d0*z^2 - 20*d0*z^3 + 15*d0*z^4 - 6*d0*z^5 + d0*z^6 + c1*c2*d2 + 5*c1*d1*z - 10*c1*d1*z^2 + 10*c1*d1*z^3 - 5*c1*d1*z^4 + c1*d1*z^5 - 4*c1*c2*d2*z + 6*c1*c2*d2*z^2 - 4*c1*c2*d2*z^3 + c1*c2*d2*z^4 - c1*c2*c3*d3 + c1*c2*c3*c4*d4 + 3*c1*c2*c3*d3*z - 3*c1*c2*c3*d3*z^2 + c1*c2*c3*d3*z^3 + c1*c2*c3*c4*d4*z^2 - c1*c2*c3*c4*c5*d5 - 2*c1*c2*c3*c4*d4*z + c1*c2*c3*c4*c5*c6*d6 + c1*c2*c3*c4*c5*d5*z)/(c1*c2 - 6*z - c1 + 5*c1*z - 10*c1*z^2 + 10*c1*z^3 - 5*c1*z^4 + c1*z^5 + 15*z^2 - 20*z^3 + 15*z^4 - 6*z^5 + z^6 - c1*c2*c3 - 4*c1*c2*z + 6*c1*c2*z^2 - 4*c1*c2*z^3 + c1*c2*z^4 + 3*c1*c2*c3*z - 3*c1*c2*c3*z^2 + c1*c2*c3*z^3 + c1*c2*c3*c4 - c1*c2*c3*c4*c5 - 2*c1*c2*c3*c4*z + c1*c2*c3*c4*z^2 + c1*c2*c3*c4*c5*c6 + c1*c2*c3*c4*c5*z + 1);
tfSev = (c1*d1 - d0 + 7*d0*z - 21*d0*z^2 + 35*d0*z^3 - 35*d0*z^4 + 21*d0*z^5 - 7*d0*z^6 + d0*z^7 - c1*c2*d2 - 6*c1*d1*z + 15*c1*d1*z^2 - 20*c1*d1*z^3 + 15*c1*d1*z^4 - 6*c1*d1*z^5 + c1*d1*z^6 + 5*c1*c2*d2*z - 10*c1*c2*d2*z^2 + 10*c1*c2*d2*z^3 - 5*c1*c2*d2*z^4 + c1*c2*d2*z^5 + c1*c2*c3*d3 - c1*c2*c3*c4*d4 - 4*c1*c2*c3*d3*z + 6*c1*c2*c3*d3*z^2 - 4*c1*c2*c3*d3*z^3 + c1*c2*c3*d3*z^4 - 3*c1*c2*c3*c4*d4*z^2 + c1*c2*c3*c4*d4*z^3 + c1*c2*c3*c4*c5*d5 + 3*c1*c2*c3*c4*d4*z - c1*c2*c3*c4*c5*c6*d6 - 2*c1*c2*c3*c4*c5*d5*z + c1*c2*c3*c4*c5*d5*z^2 + c1*c2*c3*c4*c5*c6*c7*d7 + c1*c2*c3*c4*c5*c6*d6*z)/(c1 + 7*z - c1*c2 - 6*c1*z + 15*c1*z^2 - 20*c1*z^3 + 15*c1*z^4 - 6*c1*z^5 + c1*z^6 - 21*z^2 + 35*z^3 - 35*z^4 + 21*z^5 - 7*z^6 + z^7 + c1*c2*c3 + 5*c1*c2*z - 10*c1*c2*z^2 + 10*c1*c2*z^3 - 5*c1*c2*z^4 + c1*c2*z^5 - 4*c1*c2*c3*z + 6*c1*c2*c3*z^2 - 4*c1*c2*c3*z^3 + c1*c2*c3*z^4 - c1*c2*c3*c4 + c1*c2*c3*c4*c5 + 3*c1*c2*c3*c4*z - 3*c1*c2*c3*c4*z^2 + c1*c2*c3*c4*z^3 + c1*c2*c3*c4*c5*z^2 - c1*c2*c3*c4*c5*c6 - 2*c1*c2*c3*c4*c5*z + c1*c2*c3*c4*c5*c6*c7 + c1*c2*c3*c4*c5*c6*z - 1);
tfEig = (d0 - c1*d1 - 8*d0*z + 28*d0*z^2 - 56*d0*z^3 + 70*d0*z^4 - 56*d0*z^5 + 28*d0*z^6 - 8*d0*z^7 + d0*z^8 + c1*c2*d2 + 7*c1*d1*z - 21*c1*d1*z^2 + 35*c1*d1*z^3 - 35*c1*d1*z^4 + 21*c1*d1*z^5 - 7*c1*d1*z^6 + c1*d1*z^7 - 6*c1*c2*d2*z + 15*c1*c2*d2*z^2 - 20*c1*c2*d2*z^3 + 15*c1*c2*d2*z^4 - 6*c1*c2*d2*z^5 + c1*c2*d2*z^6 - c1*c2*c3*d3 + c1*c2*c3*c4*d4 + 5*c1*c2*c3*d3*z - 10*c1*c2*c3*d3*z^2 + 10*c1*c2*c3*d3*z^3 - 5*c1*c2*c3*d3*z^4 + c1*c2*c3*d3*z^5 + 6*c1*c2*c3*c4*d4*z^2 - 4*c1*c2*c3*c4*d4*z^3 + c1*c2*c3*c4*d4*z^4 - c1*c2*c3*c4*c5*d5 - 4*c1*c2*c3*c4*d4*z + c1*c2*c3*c4*c5*c6*d6 + 3*c1*c2*c3*c4*c5*d5*z - 3*c1*c2*c3*c4*c5*d5*z^2 + c1*c2*c3*c4*c5*d5*z^3 - c1*c2*c3*c4*c5*c6*c7*d7 - 2*c1*c2*c3*c4*c5*c6*d6*z + c1*c2*c3*c4*c5*c6*d6*z^2 + c1*c2*c3*c4*c5*c6*c7*c8*d8 + c1*c2*c3*c4*c5*c6*c7*d7*z)/(c1*c2 - 8*z - c1 + 7*c1*z - 21*c1*z^2 + 35*c1*z^3 - 35*c1*z^4 + 21*c1*z^5 - 7*c1*z^6 + c1*z^7 + 28*z^2 - 56*z^3 + 70*z^4 - 56*z^5 + 28*z^6 - 8*z^7 + z^8 - c1*c2*c3 - 6*c1*c2*z + 15*c1*c2*z^2 - 20*c1*c2*z^3 + 15*c1*c2*z^4 - 6*c1*c2*z^5 + c1*c2*z^6 + 5*c1*c2*c3*z - 10*c1*c2*c3*z^2 + 10*c1*c2*c3*z^3 - 5*c1*c2*c3*z^4 + c1*c2*c3*z^5 + c1*c2*c3*c4 - c1*c2*c3*c4*c5 - 4*c1*c2*c3*c4*z + 6*c1*c2*c3*c4*z^2 - 4*c1*c2*c3*c4*z^3 + c1*c2*c3*c4*z^4 - 3*c1*c2*c3*c4*c5*z^2 + c1*c2*c3*c4*c5*z^3 + c1*c2*c3*c4*c5*c6 + 3*c1*c2*c3*c4*c5*z - c1*c2*c3*c4*c5*c6*c7 - 2*c1*c2*c3*c4*c5*c6*z + c1*c2*c3*c4*c5*c6*z^2 + c1*c2*c3*c4*c5*c6*c7*c8 + c1*c2*c3*c4*c5*c6*c7*z + 1);
tf3 = collect(simplify(tfS), z);
tf4 = collect(simplify(tfT), z);
tf5 = collect(simplify(tfF), z);
tf6 = collect(simplify(tfFive), z);
tf7 = collect(simplify(tfSix), z);
tf8 = collect(simplify(tfSev), z);
tf9 = collect(simplify(tfEig), z);
% syms b0 b1 b2 b3 b4 b5 b6 b7 b8 a1 a2 a3 a4 a5 a6 a7 a8
z1 = cplxFreq.^-1;
z2 = cplxFreq.^-2;
z3 = cplxFreq.^-3;
z4 = cplxFreq.^-4;
z5 = cplxFreq.^-5;
z6 = cplxFreq.^-6;
z7 = cplxFreq.^-7;
z8 = cplxFreq.^-8;
%% Second order
[b, a] = cheby1(2, 10, 0.6, 'low');
b0 = b(1);b1 = b(2);b2 = b(3);
a1 = a(2);a2 = a(3);
c1 = a1 + 2;
c2 = (1 + a1 + a2) / c1;
d0 = b0;
d1 = (2 * b0 + b1) / c1;
d2 = (b0 + b1 + b2) / (c1 * c2);
tfS = ((d0 - c1*d1 + c1*c2*d2)*z2 + (c1*d1 - 2*d0)*z1 + d0)./((c1*c2 - c1 + 1)*z2 + (c1 - 2)*z1 + 1);
%% Third order
[b, a] = cheby1(3, 10, 0.6, 'low');
b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);
a1 = a(2);a2 = a(3);a3 = a(4);
c1 = a1 + 3;
c2 = (3 + 2 * a1 + a2) / c1;
c3 = (1 + a1 + a2 + a3) / (c1 * c2);
d0 = b0;
d1 = (3 * b0 + b1) / c1;
d2 = (3 * b0 + 2 * b1 + b2) / (c1 * c2);
d3 = (b0 + b1 + b2 + b3) / (c1 * c2 * c3);
tfT = (d0 - 3*d0*z1 + 3*d0*z2 - d0*z3 + c1*d1*z1 - 2*c1*d1*z2 + c1*d1*z3 + c1*c2*d2*z2 - c1*c2*d2*z3 + c1*c2*c3*d3*z3)./(c1*z1 - 3*z1 - 2*c1*z2 + c1*z3 + 3*z2 - z3 + c1*c2*z2 - c1*c2*z3 + c1*c2*c3*z3 + 1);
% kDelta = [1; zeros(16383, 1)];
% out = zeros(size(kDelta));
% z1 = 0; z2 = 0; z3 = 0;
% for idx = 1 : length(kDelta)
%     x = kDelta(idx) - z1 - z2 - z3;
%     out(idx) = d0 * x + d1 * z1 + d2 * z2 + d3 * z3;
%     z3 = z3 + c3 * z2;
%     z2 = z2 + c2 * z1;
%     z1 = z1 + c1 * x;
% end
% fvtool(out)
%% Fourth order
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
tfF = (d0 - 4*d0*z1 + 6*d0*z2 - 4*d0*z3 + d0*z4 + c1*d1*z1 - 3*c1*d1*z2 + 3*c1*d1*z3 - c1*d1*z4 + c1*c2*d2*z2 - 2*c1*c2*d2*z3 + c1*c2*d2*z4 + c1*c2*c3*d3*z3 - c1*c2*c3*d3*z4 + c1*c2*c3*c4*d4*z4)./(c1*z1 - 4*z1 - 3*c1*z2 + 3*c1*z3 - c1*z4 + 6*z2 - 4*z3 + z4 + c1*c2*z2 - 2*c1*c2*z3 + c1*c2*z4 + c1*c2*c3*z3 - c1*c2*c3*z4 + c1*c2*c3*c4*z4 + 1);
% kDelta = [1; zeros(16383, 1)];
% out = zeros(size(kDelta));
% z1 = 0; z2 = 0; z3 = 0; z4 = 0;
% for idx = 1 : length(kDelta)
%     x = kDelta(idx) - z1 - z2 - z3 - z4;
%     out(idx) = d0 * x + d1 * z1 + d2 * z2 + d3 * z3 + d4 * z4;
%     z4 = z4 + c4 * z3;
%     z3 = z3 + c3 * z2;
%     z2 = z2 + c2 * z1;
%     z1 = z1 + c1 * x;
%     %% Update
%     % fc = [0.1, 0.6] + idx / length(kDelta) * 0.5;
%     % fc(fc > 0.99) = 0.99;
%     % % [b, a] = cheby1(2, 10, fc, 'bandpass');
%     % [b,a] = ellip(4,10,50,fc(1),'low');
%     % b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);b4 = b(5);
%     % a1 = a(2);a2 = a(3);a3 = a(4);a4 = a(5);
%     % c1 = a1 + 4;
%     % c2 = (6 + 3 * a1 + a2) / c1;
%     % c3 = (4 + 3 * a1 + 2 * a2 + a3) / (c1 * c2);
%     % c4 = (1 + a1 + a2 + a3 + a4) / (c1 * c2 * c3);
%     % d0 = b0;
%     % d1 = (4 * b0 + b1) / c1;
%     % d2 = (6 * b0 + 3 * b1 + b2) / (c1 * c2);
%     % d3 = (4 * b0 + 3 * b1 + 2 * b2 + b3) / (c1 * c2 * c3);
%     % d4 = (b0 + b1 + b2 + b3 + b4) / (c1 * c2 * c3 * c4);
% end
% fvtool(out)
% [h1, w1] = freqz(b,a,16384);
% [h2, w2] = freqz(out,1,16384);
% plot(w1, abs(h1));hold on;plot(w2, abs(h2));hold off;axis tight
%% Fifth order
[b,a] = ellip(5,10,50,0.5,'low');
b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);b4 = b(5);b5 = b(6);
a1 = a(2);a2 = a(3);a3 = a(4);a4 = a(5);a5 = a(6);
c1 = a1 + 5;
c2 = (10 + 4 * a1 + a2) / c1;
c3 = (10 + 6 * a1 + 3 * a2 + a3) / (c1 * c2);
c4 = (5 + 4 * a1 + 3 * a2 + 2 * a3 + a4) / (c1 * c2 * c3);
c5 = (1 + a1 + a2 + a3 + a4 + a5) / (c1 * c2 * c3 * c4);
d0 = b0;
d1 = (5 * b0 + b1) / c1;
d2 = (10 * b0 + 4 * b1 + b2) / (c1 * c2);
d3 = (10 * b0 + 6 * b1 + 3 * b2 + b3) / (c1 * c2 * c3);
d4 = (5 * b0 + 4 * b1 + 3 * b2 + 2 * b3 + b4) / (c1 * c2 * c3 * c4);
d5 = (b0 + b1 + b2 + b3 + b4 + b5) / (c1 * c2 * c3 * c4 * c5);
tfFive = (d0 - 5*d0*z1 + 10*d0*z2 - 10*d0*z3 + 5*d0*z4 - d0*z5 + c1*d1*z1 - 4*c1*d1*z2 + 6*c1*d1*z3 - 4*c1*d1*z4 + c1*d1*z5 + c1*c2*d2*z2 - 3*c1*c2*d2*z3 + 3*c1*c2*d2*z4 - c1*c2*d2*z5 + c1*c2*c3*d3*z3 - 2*c1*c2*c3*d3*z4 + c1*c2*c3*d3*z5 + c1*c2*c3*c4*d4*z4 - c1*c2*c3*c4*d4*z5 + c1*c2*c3*c4*c5*d5*z5)./(c1*z1 - 5*z1 - 4*c1*z2 + 6*c1*z3 - 4*c1*z4 + c1*z5 + 10*z2 - 10*z3 + 5*z4 - z5 + c1*c2*z2 - 3*c1*c2*z3 + 3*c1*c2*z4 - c1*c2*z5 + c1*c2*c3*z3 - 2*c1*c2*c3*z4 + c1*c2*c3*z5 + c1*c2*c3*c4*z4 - c1*c2*c3*c4*z5 + c1*c2*c3*c4*c5*z5 + 1);
%% Sixth order
[b,a] = ellip(6,10,50,0.5,'low');
[b, a] = ellip(3, 10, 50, [0.3, 0.6], 'bandpass');
b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);b4 = b(5);b5 = b(6);b6 = b(7);
a1 = a(2);a2 = a(3);a3 = a(4);a4 = a(5);a5 = a(6);a6 = a(7);
c1 = a1 + 6;
c2 = (15 + 5 * a1 + a2) / c1;
c3 = (20 + 10 * a1 + 4 * a2 + a3) / (c1 * c2);
c4 = (15 + 10 * a1 + 6 * a2 + 3 * a3 + a4) / (c1 * c2 * c3);
c5 = (6 + 5 * a1 + 4 * a2 + 3 * a3 + 2 * a4 + a5) / (c1 * c2 * c3 * c4);
c6 = (1 + a1 + a2 + a3 + a4 + a5 + a6) / (c1 * c2 * c3 * c4 * c5);
d0 = b0;
d1 = (6 * b0 + b1) / c1;
d2 = (15 * b0 + 5 * b1 + b2) / (c1 * c2);
d3 = (20 * b0 + 10 * b1 + 4 * b2 + b3) / (c1 * c2 * c3);
d4 = (15 * b0 + 10 * b1 + 6 * b2 + 3 * b3 + b4) / (c1 * c2 * c3 * c4);
d5 = (6 * b0 + 5 * b1 + 4 * b2 + 3 * b3 + 2 * b4 + b5) / (c1 * c2 * c3 * c4 * c5);
d6 = (b0 + b1 + b2 + b3 + b4 + b5 + b6) / (c1 * c2 * c3 * c4 * c5 * c6);
tfSix = (d0 - 6*d0*z1 + 15*d0*z2 - 20*d0*z3 + 15*d0*z4 - 6*d0*z5 + d0*z6 + c1*d1*z1 - 5*c1*d1*z2 + 10*c1*d1*z3 - 10*c1*d1*z4 + 5*c1*d1*z5 - c1*d1*z6 + c1*c2*d2*z2 - 4*c1*c2*d2*z3 + 6*c1*c2*d2*z4 - 4*c1*c2*d2*z5 + c1*c2*d2*z6 + c1*c2*c3*d3*z3 - 3*c1*c2*c3*d3*z4 + 3*c1*c2*c3*d3*z5 - c1*c2*c3*d3*z6 + c1*c2*c3*c4*d4*z4 - 2*c1*c2*c3*c4*d4*z5 + c1*c2*c3*c4*d4*z6 + c1*c2*c3*c4*c5*d5*z5 - c1*c2*c3*c4*c5*d5*z6 + c1*c2*c3*c4*c5*c6*d6*z6)./(c1*z1 - 6*z1 - 5*c1*z2 + 10*c1*z3 - 10*c1*z4 + 5*c1*z5 - c1*z6 + 15*z2 - 20*z3 + 15*z4 - 6*z5 + z6 + c1*c2*z2 - 4*c1*c2*z3 + 6*c1*c2*z4 - 4*c1*c2*z5 + c1*c2*z6 + c1*c2*c3*z3 - 3*c1*c2*c3*z4 + 3*c1*c2*c3*z5 - c1*c2*c3*z6 + c1*c2*c3*c4*z4 - 2*c1*c2*c3*c4*z5 + c1*c2*c3*c4*z6 + c1*c2*c3*c4*c5*z5 - c1*c2*c3*c4*c5*z6 + c1*c2*c3*c4*c5*c6*z6 + 1);
% plot(abs(tfSix));axis tight
%% Seventh order
[b,a] = ellip(7,10,50,0.5,'low');
b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);b4 = b(5);b5 = b(6);b6 = b(7);b7 = b(8);
a1 = a(2);a2 = a(3);a3 = a(4);a4 = a(5);a5 = a(6);a6 = a(7);a7 = a(8);
c1 = a1 + 7;
c2 = (21 + 6 * a1 + a2) / c1;
c3 = (35 + 15 * a1 + 5 * a2 + a3) / (c1 * c2);
c4 = (35 + 20 * a1 + 10 * a2 + 4 * a3 + a4) / (c1 * c2 * c3);
c5 = (21 + 15 * a1 + 10 * a2 + 6 * a3 + 3 * a4 + a5) / (c1 * c2 * c3 * c4);
c6 = (7 + 6 * a1 + 5 * a2 + 4 * a3 + 3 * a4 + 2 * a5 + a6) / (c1 * c2 * c3 * c4 * c5);
c7 = (1 + a1 + a2 + a3 + a4 + a5 + a6 + a7) / (c1 * c2 * c3 * c4 * c5 * c6);
d0 = b0;
d1 = (7 * b0 + b1) / c1;
d2 = (21 * b0 + 6 * b1 + b2) / (c1 * c2);
d3 = (35 * b0 + 15 * b1 + 5 * b2 + b3) / (c1 * c2 * c3);
d4 = (35 * b0 + 20 * b1 + 10 * b2 + 4 * b3 + b4) / (c1 * c2 * c3 * c4);
d5 = (21 * b0 + 15 * b1 + 10 * b2 + 6 * b3 + 3 * b4 + b5) / (c1 * c2 * c3 * c4 * c5);
d6 = (7 * b0 + 6 * b1 + 5 * b2 + 4 * b3 + 3 * b4 + 2 * b5 + b6) / (c1 * c2 * c3 * c4 * c5 * c6);
d7 = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7) / (c1 * c2 * c3 * c4 * c5 * c6 * c7);
tfSev = (d0 - 7*d0*z1 + 21*d0*z2 - 35*d0*z3 + 35*d0*z4 - 21*d0*z5 + 7*d0*z6 - d0*z7 + c1*d1*z1 - 6*c1*d1*z2 + 15*c1*d1*z3 - 20*c1*d1*z4 + 15*c1*d1*z5 - 6*c1*d1*z6 + c1*d1*z7 + c1*c2*d2*z2 - 5*c1*c2*d2*z3 + 10*c1*c2*d2*z4 - 10*c1*c2*d2*z5 + 5*c1*c2*d2*z6 - c1*c2*d2*z7 + c1*c2*c3*d3*z3 - 4*c1*c2*c3*d3*z4 + 6*c1*c2*c3*d3*z5 - 4*c1*c2*c3*d3*z6 + c1*c2*c3*d3*z7 + c1*c2*c3*c4*d4*z4 - 3*c1*c2*c3*c4*d4*z5 + 3*c1*c2*c3*c4*d4*z6 - c1*c2*c3*c4*d4*z7 + c1*c2*c3*c4*c5*d5*z5 - 2*c1*c2*c3*c4*c5*d5*z6 + c1*c2*c3*c4*c5*d5*z7 + c1*c2*c3*c4*c5*c6*d6*z6 - c1*c2*c3*c4*c5*c6*d6*z7 + c1*c2*c3*c4*c5*c6*c7*d7*z7)./(c1*z1 - 7*z1 - 6*c1*z2 + 15*c1*z3 - 20*c1*z4 + 15*c1*z5 - 6*c1*z6 + c1*z7 + 21*z2 - 35*z3 + 35*z4 - 21*z5 + 7*z6 - z7 + c1*c2*z2 - 5*c1*c2*z3 + 10*c1*c2*z4 - 10*c1*c2*z5 + 5*c1*c2*z6 - c1*c2*z7 + c1*c2*c3*z3 - 4*c1*c2*c3*z4 + 6*c1*c2*c3*z5 - 4*c1*c2*c3*z6 + c1*c2*c3*z7 + c1*c2*c3*c4*z4 - 3*c1*c2*c3*c4*z5 + 3*c1*c2*c3*c4*z6 - c1*c2*c3*c4*z7 + c1*c2*c3*c4*c5*z5 - 2*c1*c2*c3*c4*c5*z6 + c1*c2*c3*c4*c5*z7 + c1*c2*c3*c4*c5*c6*z6 - c1*c2*c3*c4*c5*c6*z7 + c1*c2*c3*c4*c5*c6*c7*z7 + 1);
%% Eighth order
[b,a] = ellip(8,10,50,0.5,'low');
[b, a] = ellip(4, 10, 50, [0.3, 0.6], 'bandpass');
b0 = b(1);b1 = b(2);b2 = b(3);b3 = b(4);b4 = b(5);b5 = b(6);b6 = b(7);b7 = b(8);b8 = b(9);
a1 = a(2);a2 = a(3);a3 = a(4);a4 = a(5);a5 = a(6);a6 = a(7);a7 = a(8);a8 = a(9);
c1 = a1 + 8;
c2 = (28 + 7 * a1 + a2) / c1;
c3 = (56 + 21 * a1 + 6 * a2 + a3) / (c1 * c2);
c4 = (70 + 35 * a1 + 15 * a2 + 5 * a3 + a4) / (c1 * c2 * c3);
c5 = (56 + 35 * a1 + 20 * a2 + 10 * a3 + 4 * a4 + a5) / (c1 * c2 * c3 * c4);
c6 = (28 + 21 * a1 + 15 * a2 + 10 * a3 + 6 * a4 + 3 * a5 + a6) / (c1 * c2 * c3 * c4 * c5);
c7 = (8 + 7 * a1 + 6 * a2 + 5 * a3 + 4 * a4 + 3 * a5 + 2 * a6 + a7) / (c1 * c2 * c3 * c4 * c5 * c6);
c8 = (1 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8) / (c1 * c2 * c3 * c4 * c5 * c6 * c7);
d0 = b0;
d1 = (8 * b0 + b1) / c1;
d2 = (28 * b0 + 7 * b1 + b2) / (c1 * c2);
d3 = (56 * b0 + 21 * b1 + 6 * b2 + b3) / (c1 * c2 * c3);
d4 = (70 * b0 + 35 * b1 + 15 * b2 + 5 * b3 + b4) / (c1 * c2 * c3 * c4);
d5 = (56 * b0 + 35 * b1 + 20 * b2 + 10 * b3 + 4 * b4 + b5) / (c1 * c2 * c3 * c4 * c5);
d6 = (28 * b0 + 21 * b1 + 15 * b2 + 10 * b3 + 6 * b4 + 3 * b5 + b6) / (c1 * c2 * c3 * c4 * c5 * c6);
d7 = (8 * b0 + 7 * b1 + 6 * b2 + 5 * b3 + 4 * b4 + 3 * b5 + 2 * b6 + b7) / (c1 * c2 * c3 * c4 * c5 * c6 * c7);
d8 = (b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8) / (c1 * c2 * c3 * c4 * c5 * c6 * c7 * c8);
tfEig = (d0 - 8*d0*z1 + 28*d0*z2 - 56*d0*z3 + 70*d0*z4 - 56*d0*z5 + 28*d0*z6 - 8*d0*z7 + d0*z8 + c1*d1*z1 - 7*c1*d1*z2 + 21*c1*d1*z3 - 35*c1*d1*z4 + 35*c1*d1*z5 - 21*c1*d1*z6 + 7*c1*d1*z7 - c1*d1*z8 + c1*c2*d2*z2 - 6*c1*c2*d2*z3 + 15*c1*c2*d2*z4 - 20*c1*c2*d2*z5 + 15*c1*c2*d2*z6 - 6*c1*c2*d2*z7 + c1*c2*d2*z8 + c1*c2*c3*d3*z3 - 5*c1*c2*c3*d3*z4 + 10*c1*c2*c3*d3*z5 - 10*c1*c2*c3*d3*z6 + 5*c1*c2*c3*d3*z7 - c1*c2*c3*d3*z8 + c1*c2*c3*c4*d4*z4 - 4*c1*c2*c3*c4*d4*z5 + 6*c1*c2*c3*c4*d4*z6 - 4*c1*c2*c3*c4*d4*z7 + c1*c2*c3*c4*d4*z8 + c1*c2*c3*c4*c5*d5*z5 - 3*c1*c2*c3*c4*c5*d5*z6 + 3*c1*c2*c3*c4*c5*d5*z7 - c1*c2*c3*c4*c5*d5*z8 + c1*c2*c3*c4*c5*c6*d6*z6 - 2*c1*c2*c3*c4*c5*c6*d6*z7 + c1*c2*c3*c4*c5*c6*d6*z8 + c1*c2*c3*c4*c5*c6*c7*d7*z7 - c1*c2*c3*c4*c5*c6*c7*d7*z8 + c1*c2*c3*c4*c5*c6*c7*c8*d8*z8)./(c1*z1 - 8*z1 - 7*c1*z2 + 21*c1*z3 - 35*c1*z4 + 35*c1*z5 - 21*c1*z6 + 7*c1*z7 - c1*z8 + 28*z2 - 56*z3 + 70*z4 - 56*z5 + 28*z6 - 8*z7 + z8 + c1*c2*z2 - 6*c1*c2*z3 + 15*c1*c2*z4 - 20*c1*c2*z5 + 15*c1*c2*z6 - 6*c1*c2*z7 + c1*c2*z8 + c1*c2*c3*z3 - 5*c1*c2*c3*z4 + 10*c1*c2*c3*z5 - 10*c1*c2*c3*z6 + 5*c1*c2*c3*z7 - c1*c2*c3*z8 + c1*c2*c3*c4*z4 - 4*c1*c2*c3*c4*z5 + 6*c1*c2*c3*c4*z6 - 4*c1*c2*c3*c4*z7 + c1*c2*c3*c4*z8 + c1*c2*c3*c4*c5*z5 - 3*c1*c2*c3*c4*c5*z6 + 3*c1*c2*c3*c4*c5*z7 - c1*c2*c3*c4*c5*z8 + c1*c2*c3*c4*c5*c6*z6 - 2*c1*c2*c3*c4*c5*c6*z7 + c1*c2*c3*c4*c5*c6*z8 + c1*c2*c3*c4*c5*c6*c7*z7 - c1*c2*c3*c4*c5*c6*c7*z8 + c1*c2*c3*c4*c5*c6*c7*c8*z8 + 1);
