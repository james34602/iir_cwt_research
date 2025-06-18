function Deriche4nd_svf(x, sigma)
if nargin == 0
    sigma = 0.1;
    sigLen = round(sigma * 10);
    if mod(sigLen, 2) == 1
        sigLen = sigLen + 1;
    end
    if sigLen < 16
        sigLen = 16;
    end
    x = zeros(sigLen, 1);
    x(sigLen / 2 + 1) = 1;
end
numSigmas = length(x);
numSigmas = 1025;rng(1);x = randn(numSigmas,1);
sigma = linspace(4, 30, numSigmas) .^ 4;
sigma = sigma ./ max(sigma) * 30 + 5;
% sigma = rand(numSigmas, 1) * 30;[b,a]=butter(4, 0.04);sigma = filtfilt(b,a,sigma);
% sigma(:) = 0.1;
[bfwd, bbwd, a, Bmfwd, Bmbwd, Am, sosfwd, sosbwd] = InitDeriche(sigma(1), numSigmas);
bfwd_ = zeros([size(bfwd), numSigmas]);
bbwd_ = zeros([size(bbwd), numSigmas]);
a_ = zeros([size(a), numSigmas]);
Bmfwd_ = zeros([size(Bmfwd), numSigmas]);
Bmbwd_ = zeros([size(Bmbwd), numSigmas]);
Am_ = zeros([size(Am), numSigmas]);
sosfwd_ = zeros([size(sosfwd), numSigmas]);
sosbwd_ = zeros([size(sosbwd), numSigmas]);
divareaSum = zeros(numSigmas, 1);
for i = 1 : numSigmas
    [bfwd_(:, :, i), bbwd_(:, :, i), a_(:, :, i), Bmfwd_(:, :, i), Bmbwd_(:, :, i), Am_(:, :, i), sosfwd_(:, :, i), sosbwd_(:, :, i), divareaSum(i)] = InitDeriche(sigma(i), numSigmas);
end
%% LTV
modrespfwd_=smpparsvffilt2(Bmfwd_, Am_, x);
Bmbwd = cat(3, Bmbwd_(:, :, 1), Bmbwd_);
Am = cat(3, Am_(:, :, 1), Am_);
Bmbwd = flip(Bmbwd, 3);
Am = flip(Am, 3);
modrespbwd_=flipud(smpparsvffilt2(Bmbwd, Am, [0; flipud(x)]));
modrespbwd_(1) = [];
modresp_ = (modrespfwd_ + modrespbwd_) .* divareaSum;
y = smoothSpectrum2(x, sigma);
bk2 = smpparsvffiltbidirectional2(Bmfwd_, Bmbwd_, Am_, x) .* divareaSum;
bk1 = smpparsvffiltbidirectional1(Bmfwd_, Bmbwd_, Am_, x) .* divareaSum;
plot(y);
hold on
plot(bk2)
hold off;
axis tight
return
%% LTV
resfwd_=smpsossvffilt2(sosfwd_,x);
sosbwd_ = cat(3, sosbwd_(:, :, 1), sosbwd_);
sosbwd_ = flip(sosbwd_, 3);
resbwd_=flipud(smpsossvffilt2(sosbwd_,[0; flipud(x)]));
resbwd_(1) = [];
res_ = (resfwd_ + resbwd_) .* divareaSum;
plot(modresp_);
hold on
plot(res_);
hold off
axis tight
end
function y = smpparsvffiltbidirectional2(b1, b2, a, x)
c1 = a(2, :, :) + 2;
c2 = (1 + a(2, :, :) + a(3, :, :)) ./ c1;
d0fwd = b1(1, :, :);
d1fwd = (2 * b1(1, :, :) + b1(2, :, :)) ./ c1;
d2fwd = (b1(1, :, :) + b1(2, :, :)) ./ (c1 .* c2);
d0bwd = b2(1, :, :);
d1bwd = (2 * b2(1, :, :) + b2(2, :, :)) ./ c1;
d2bwd = (b2(1, :, :) + b2(2, :, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y1 = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0fwd(1, 1, i) * Y + d1fwd(1, 1, i) * z1_A(1) + d2fwd(1, 1, i) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0fwd(1, 2, i) * Y2 + d1fwd(1, 2, i) * z1_A(2) + d2fwd(1, 2, i) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i) * Y2;
    y1(i) = Yi1 + Yi2;
end
z1_A = [0, 0];
z2_A = [0, 0];
y2 = zeros(size(x, 1), 1);
for i = length(x) : -1 : 2
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0bwd(1, 1, i - 1) * Y + d1bwd(1, 1, i - 1) * z1_A(1) + d2bwd(1, 1, i - 1) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i - 1) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i - 1) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0bwd(1, 2, i - 1) * Y2 + d1bwd(1, 2, i - 1) * z1_A(2) + d2bwd(1, 2, i - 1) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i - 1) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i - 1) * Y2;
    y2(i - 1) = Yi1 + Yi2;
end
y = y1 + y2;
end
function y1 = smpparsvffiltbidirectional1(b1, b2, a, x)
c1 = a(2, :, :) + 2;
c2 = (1 + a(2, :, :) + a(3, :, :)) ./ c1;
d0fwd = b1(1, :, :);
d1fwd = (2 * b1(1, :, :) + b1(2, :, :)) ./ c1;
d2fwd = (b1(1, :, :) + b1(2, :, :)) ./ (c1 .* c2);
d0bwd = b2(1, :, :);
d1bwd = (2 * b2(1, :, :) + b2(2, :, :)) ./ c1;
d2bwd = (b2(1, :, :) + b2(2, :, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y1 = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0fwd(1, 1, i) * Y + d1fwd(1, 1, i) * z1_A(1) + d2fwd(1, 1, i) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0fwd(1, 2, i) * Y2 + d1fwd(1, 2, i) * z1_A(2) + d2fwd(1, 2, i) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i) * Y2;
    y1(i) = Yi1 + Yi2;
end
z1_A = [0, 0];
z2_A = [0, 0];
y2 = zeros(size(x, 1) - 1, 1);
for i = length(x) : -1 : 2
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0bwd(1, 1, i) * Y + d1bwd(1, 1, i) * z1_A(1) + d2bwd(1, 1, i) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0bwd(1, 2, i) * Y2 + d1bwd(1, 2, i) * z1_A(2) + d2bwd(1, 2, i) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i) * Y2;
    y2(i - 1) = Yi1 + Yi2;
end
y1(1 : (size(x, 1) - 1)) = y1(1 : (size(x, 1) - 1)) + y2;
end
function x_oct = smoothSpectrum2(X,sigma)
f = (0:(length(X) - 1))';
x_oct = zeros(size(X));
for i = 1 : length(sigma)
    g = exp(-( ( (f - f(i) ) .^ 2 ) / (2 * (sigma(i) ^ 2) ) ) ); % Gaussian
    g = g ./ sum(g);
    x_oct(i) = sum(g.*X); % calculate smoothed spectral coefficient
end
end
function y = smpparsvffilt2(b, a, x)
c1 = a(2, :, :) + 2;
c2 = (1 + a(2, :, :) + a(3, :, :)) ./ c1;
d0 = b(1, :, :);
d1 = (2 * b(1, :, :) + b(2, :, :)) ./ c1;
d2 = (b(1, :, :) + b(2, :, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0(1, 1, i) * Y + d1(1, 1, i) * z1_A(1) + d2(1, 1, i) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0(1, 2, i) * Y2 + d1(1, 2, i) * z1_A(2) + d2(1, 2, i) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i) * Y2;
    y(i) = Yi1 + Yi2;
end
end
function y = smpsossvffilt2(sos, x)
b = sos(:, 1 : 3, :);
a = sos(:, 3 + 1 : 6, :);
b = permute(b, [2, 1, 3]);
a = permute(a, [2, 1, 3]);
c1 = a(2, :, :) + 2;
c2 = (1 + a(2, :, :) + a(3, :, :)) ./ c1;
d0 = b(1, :, :);
d1 = (2 * b(1, :, :) + b(2, :, :)) ./ c1;
d2 = (b(1, :, :) + b(2, :, :) + b(3, :, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0(1, 1, i) * Y + d1(1, 1, i) * z1_A(1) + d2(1, 1, i) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i) * Y;
    Y2 = Yi1 - z1_A(2) - z2_A(2);
    y(i) = d0(1, 2, i) * Y2 + d1(1, 2, i) * z1_A(2) + d2(1, 2, i) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i) * Y2;
end
end
function y = smpparsvffilt(b, a, x)
c1 = a(2, :) + 2;
c2 = (1 + a(2, :) + a(3, :)) ./ c1;
d0 = b(1, :);
d1 = (2 * b(1, :) + b(2, :)) ./ c1;
d2 = (b(1, :) + b(2, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0(1) * Y + d1(1) * z1_A(1) + d2(1) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0(2) * Y2 + d1(2) * z1_A(2) + d2(2) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(2) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(2) * Y2;
    y(i) = Yi1 + Yi2;
end
end
function y = smpsossvffilt(sos, x)
b = sos(:, 1 : 3)';
a = sos(:, 3 + 1 : 6)';
c1 = a(2, :) + 2;
c2 = (1 + a(2, :) + a(3, :)) ./ c1;
d0 = b(1, :);
d1 = (2 * b(1, :) + b(2, :)) ./ c1;
d2 = (b(1, :) + b(2, :) + b(3, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0(1) * Y + d1(1) * z1_A(1) + d2(1) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1) * Y;
    Y2 = Yi1 - z1_A(2) - z2_A(2);
    y(i) = d0(2) * Y2 + d1(2) * z1_A(2) + d2(2) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(2) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(2) * Y2;
end
end