halfLen = 4097;
halfWndLen = halfLen - 1;
digw = linspace(0, pi - pi / halfWndLen, halfWndLen);
digw(halfWndLen) = pi - pi / halfWndLen;
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
[b, a] = cheby1(1, 10, [0.6, 0.8], 'bandpass');
% [b, a] = cheby1(2, 10, 0.6, 'low');
% [b, a] = cheby1(2, 10, 0.6, 'high');
b0 = b(1);b1 = b(2);b2 = b(3);
a1 = a(2);a2 = a(3);
z1 = cplxFreq.^-1;
z2 = cplxFreq.^-2;
% %% Arbitrary
c1 = a1 + 2;
c2 = (1 + a1 + a2) / c1;
d0 = b0;
d1 = (2 * b0 + b1) / c1;
d2 = (b0 + b1 + b2) / (c1 * c2);
rec_a1 = c1 - 2;
rec_a2 = c1*c2 - c1 + 1;
rec_b1 = d1 * c1 - 2 * d0;
rec_b2 = d0 - c1*d1 + c1*c2*d2;
%% Text book filter
% F = 1/24;
% Q = sqrt(2)/2;
% w = 2 + tan(pi + F);
% a = w/Q;
% b = w * w;
% c1 = (a + b)/(1 + a/2 + b/4);
% c2 = b/(a + b);
% %% High pass
% d0 = 1 - c1/2 + c1 * c2/4;
% d1 = 0;
% d2 = 0;
% %% Band pass
% d0 = (1 - c2) * c1/2;
% d1 = 1 - c2;
% d2 = 0;
% %% Low pass
% d0 = c1 * c2/4;
% d1 = c2;
% d2 = 1;
% %% Evaluate transfer
syms z a1 a2 b0 b1 b2 x y hz
z1 = z^-1;
z2 = z^-2;
c1 = a1 + 2;
c2 = (1 + a1 + a2) / c1;
d0 = b0;
d1 = (2 * b0 + b1) / c1;
d2 = (b0 + b1 + b2) / (c1 * c2);
syms c1 c2 d0 d1 d2
a1 = c1 - 2;
a2 = c1 * c2 - c1 + 1;
b0 = d0;
b1 = d1 * c1 - 2 * d0;
b2 = d0 - c1 * d1 + c1 * c2 * d2;
Hz = (b0 + b1 * z1 + b2 * z2) / (1 + a1 * z1 + a2 * z2);
term = c1*z1 - 2*z1 - c1*z2 + z2 + c1*c2*z2 + 1;
branch1 = (z2 - 2*z1 + 1) ./ term;
branch2 = (c1 * (z1 - z2)) ./ term;
branch3 = (c1*c2 * z2) ./ term;
tf = d0 * branch1 + d1 * branch2 + d2 * branch3;
tf2 = (d0 - 2*d0*z1 + d0*z2 + c1*d1*z1 - c1*d1*z2 + c1*c2*d2*z2) ./ term;
%% Evaluate transfer for second order Gaussian approximation
c1 = a1 + 2;
c2 = (1 + a1 + a2) / c1;
d0 = b0;
tf2 = (d0 - 2*d0*z1 + d0*z2 + 2*c1*c2*z1 + c1*c2*z2 - 2*c1*c2*z2) ./ term;