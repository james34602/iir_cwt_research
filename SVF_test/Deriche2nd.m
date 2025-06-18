function [y, b1, a1, b2, a2] = Deriche2nd(x,fsigma)
npixels = length(x);

alpha = 1.695 / fsigma;
e =  exp(-alpha);
e2 = exp(-2 * alpha);
norm = (1 - e) * (1 - e) / (1 + 2 * alpha * e - e2) * 0.95;

p_filter(1) = -e2;
p_filter(2) = 2 * e;
p_filter(3) = norm;
p_filter(4) = norm * (alpha - 1) * e;
p_filter(5) = 0;

n_filter(1) = -e2;
n_filter(2) = 2 * e;
n_filter(3) = 0;
n_filter(4) = norm * (alpha + 1) * e;
n_filter(5) = -norm * e2;

%
%--- Now let's apply the Deriche equation (15)
%
%    We run left to right accross the line of pixels
%

d2 = p_filter(1);
d1 = p_filter(2);

n0 = p_filter(3);
n1 = p_filter(4);
n2 = p_filter(5);   % Note this is always == 0

in1 = x(1);
in2 = x(1);

out1 = (n2 + n1 + n0)*in1 / (1.0-d1-d2);
out2 = (n2 + n1 + n0)*in1 / (1.0-d1-d2);

y1 = zeros(1, npixels);
for i=1:npixels
    in0  = x(i);
    out0 = n2*in2 + n1*in1 + n0*in0 + d1*out1 + d2*out2;
    in2  = in1;
    in1  = in0;
    out2 = out1;
    out1 = out0;
    y1(i) = out0;
end
b1 = [n0, n1, n2];
a1 = [1, -d1, -d2];
% fvtool([n0, n1, n2], [1, -d1, -d2])
% We run right to left accross the line of pixels
d2 = n_filter(1);
d1 = n_filter(2);
n0 = n_filter(3);   % Always == 0
n1 = n_filter(4);
n2 = n_filter(5);

in1 = x(npixels);
in2 = x(npixels);

out1 = (n2 + n1 + n0)*in1/(1.0-d1-d2);
out2 = (n2 + n1 + n0)*in1/(1.0-d1-d2);

y2 = zeros(1, npixels);
for i=npixels: - 1 : 1
    in0  = x(i);
    out0 = n2*in2 + n1*in1 + n0*in0 + d1*out1 + d2*out2;
    in2  = in1;
    in1  = in0;
    out2 = out1;
    out1 = out0;
    y2(i) = out0;
end
b2 = [n0, n1, n2];
a2 = [1, -d1, -d2];
% fvtool([n0, n1, n2], [1, -d1, -d2])
% The final result is the summation of the vectors produced by equations 15
% and 16 of Deriche's paper.
y = y1 + y2;
end