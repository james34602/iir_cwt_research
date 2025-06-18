addpath('../')
[b, a] = butter(1, [0.1, 0.2], 'stop');
b0=b(1);b1=b(2);b2=b(3);a2=a(2);a3=a(3);
sigLen = 3001;
centre = 1500;
kDelta = zeros(sigLen, 1);
kDelta(centre + 1) = 1;
%% TF2SVF
c1 = a2 + 2;
c2 = (1 + a2 + a3) / c1;
d0 = b0;
d1 = (2 * b0 + b1) / c1;
d2 = (b0 + b1 + b2) / (c1 * c2);
%% SVF2TF
a2 = c1 - 2;
a3 = (c2 * c1) - 1 - a2;
b0 = d0;
b1 = (d1 * c1) - 2 * d0;
b2 = (d2 * c1 * c2) - b0 - b1;
%% Vanilla SVF debugging
z1_A = 0.0;
z2_A = 0.0;
oo = zeros(size(kDelta, 1), 1);
c1 = a(2) + 2;
c2 = (1 + a(2) + a(3)) / c1;
d1 = (2 * b0 + b1) / c1;
d2 = (b0 + b1 + b2) / (c1 * c2);
for idx = 1 : size(kDelta, 1)
    x = kDelta(idx);
    y = x - z1_A - z2_A;
    oo(idx) = d0 * y + d1 * z1_A + d2 * z2_A;
    z2_A = z2_A + c2 * z1_A;
    z1_A = z1_A + c1 * y;
end
% Initialize output y and previous values
y = zeros(size(kDelta, 1), 1);
for n = 1 : size(kDelta, 1)
    if n == 1
        y(n) = d0*0 - 2*d0*0 + d0*kDelta(n) + c1*d1*0 - c1*d1*0 + c1*c2*d2*0 - c1*c2*0 + c1*0 - c1*0 + 2*0 - 0;
    elseif n == 2
        y(n) = d0*0 - 2*d0*kDelta(n - 1) + d0*kDelta(n) + c1*d1*kDelta(n - 1) - c1*d1*0 + c1*c2*d2*0 - c1*c2*0 + c1*0 - c1*y(n - 1) + 2*y(n - 1) - 0;
    else
        y(n) = d0*kDelta(n - 2) - 2*d0*kDelta(n - 1) + d0*kDelta(n) + c1*d1*kDelta(n - 1) - c1*d1*kDelta(n - 2) + c1*c2*d2*kDelta(n - 2) - c1*c2*y(n - 2) + c1*y(n - 2) - c1*y(n - 1) + 2*y(n - 1) - y(n - 2);
    end
end
plot(oo);hold on;plot(y);hold off