addpath('../')
addpath('../multResByIIRInFreq_FullyFreqDomainCorrMatrix')
sigLen = 3001;
centre = 1500;
kDelta = zeros(sigLen, 1);
kDelta(centre + 1) = 1;
%%
rng(3)
nSigmas = 6;
sigmas = [100, 200, 150, 20, 40, 3];
[b, a, c1, c2] = gauss_precompute(sigmas);
%%
sel = 2;
x = 0:1:(sigLen-1);
standardGauss = gaussmf(x, [sigmas(sel), centre]);
standardGauss = standardGauss ./ sum(standardGauss);
vYSignal = filter(b(sel), a(sel, :), kDelta);
vYSignal = filter(b(sel), a(sel, :), vYSignal(end:-1:1));
vYSignal = vYSignal(end:-1:1);
%%
cb = zeros(size(kDelta, 1), 1);cc1 = zeros(size(kDelta, 1), 1);cc2 = zeros(size(kDelta, 1), 1);
cb(:) = b(sel);
cc1(:) = c1(sel);
cc2(:) = c2(sel);
tmp = zeros(size(kDelta, 1), 1, 'like', kDelta);
% q_fft_frame = ltv1Slice(kDelta, tmp, cb, [], cc1, cc2);
[q_fft_frame2, b1, a1, b2, a2] = Deriche2nd(kDelta, sigmas(sel));
[~, ~, ~, ~, ~, ~, ~, ~, divareaSum, q_fft_frame3] = InitDeriche(sigmas(sel), size(kDelta, 1));
%%
figure(1)
plot(standardGauss)
hold on;
plot(vYSignal);
% plot(q_fft_frame);
plot(q_fft_frame2);
plot(q_fft_frame3 * divareaSum);
hold off;
axis tight;
legend('True Gaussian', 'FwdBwd 2nd', 'Deriche 2nd', 'Deriche 4nd');
title('Approximating a Gaussian')
% return
%% Gaussian sigma to time constant
% fc3db = sqrt(log(2)) ./ (2 * sigmas * pi);
% alpha = cos(fc3db) - 1 + sqrt(cos(fc3db).^2 - 4*cos(fc3db) + 3);
% time = (1 - alpha) ./ alpha;
%% 1/N octave spectrum smoothing
clear
load handel.mat
fs = 48000;
% y = zeros(2048, 1);
% y(2109) = 1;
y = y(1 : 4096);
fftLen = length(y);
if mod(fftLen,2)==0 % keep only meaningful frequencies
    halfLen = (fftLen/2)+1;
else
    halfLen = (fftLen+1)/2;
end
Y = fft(fftshift(y .* hann(fftLen)));
Y = abs(Y);
% frequency bins grid (linear in this case) - pre and pos padding is added
% poles of the IIR LTV Q FFT transform for the parameters above
f = (0:1:fftLen/2)*fs/fftLen;
f2 = f;
f2(halfLen+1:fftLen) = conj(f2(halfLen-1:-1:2));
% number of points of pre and post padding used to set initial conditions
thetas = 0:(fftLen/2);
thetas(1) = eps;
thetas(halfLen+1:fftLen) = conj(thetas(halfLen-1:-1:2));
% frequency bins grid (linear in this case) - pre and pos padding is added
Q = 3.5;
Noct = estimate1NOctFromQ(Q + (0.9/3),fftLen,3);
% thetas = [thetas(end-prepad+1:end), thetas(1 : halfLen + pospad - 1)];
%% Smoothing by time varying filter
% Noct = 8;
sigmas = (thetas ./ fftLen) ./ Noct / pi * fftLen;
[b, ~, c1, c2] = gauss_precompute(sigmas);
Z1 = smoothSpectrum2(Y, f2', Noct);
[~, ~, ~, Bmfwd, Bmbwd, Am, ~, ~, ~] = InitDeriche(sigmas(1), length(thetas));
Bmfwd = zeros([size(Bmfwd), length(thetas)]);
Bmbwd = zeros([size(Bmbwd), length(thetas)]);
Am = zeros([size(Am), length(thetas)]);
divareaSum = zeros(length(thetas), 1);
pts = 0;
for i = 1 : length(thetas)
    [~, ~, ~, Bmfwd(:, :, i), Bmbwd(:, :, i), Am(:, :, i), ~, ~, divareaSum(i)] = InitDeriche(sigmas(i), length(thetas));
end
%% First order time varying filtering
fgt_facM = 0.3 / Noct;
poles = exp(-1 ./ (fgt_facM .* thetas'));
tmp1 = zeros(length(Y), 1);
firstOrderOut = zeros(length(Y), 1);
% Forward filtering
z = 0;
for n=1:length(Y)
    o = Y(n) + z;
    z = poles(n) * o;
    tmp1(n) = o;
end
% Reverse filtering
z = 0;
for n=length(Y):-1:1
    o = tmp1(n) + z;
    z = poles(n) * o;
    firstOrderOut(n) = o;
end
firstOrderOut = firstOrderOut .* (1 - poles).^2;
%% Second order time varying filtering
tmp2 = zeros(length(Y), 1);
secondOrderOut = zeros(length(Y), 1);
% Forward filtering
z1_A = 0; z2_A = 0;
for a = 1:length(Y)
    x = Y(a);
    y = x - z1_A - z2_A;
    st = c2(a) * z1_A;
    tmp2(a) = b(a) * y + 2 * st + z2_A;
    z2_A = z2_A + st;
    z1_A = z1_A + c1(a) * y;
end
% Reverse filtering
z1_A = 0; z2_A = 0;
for a = length(Y):-1:1
    x = tmp2(a);
    y = x - z1_A - z2_A;
    st = c2(a) * z1_A;
    secondOrderOut(a) = b(a) * y + 2 * st + z2_A;
    z2_A = z2_A + st;
    z1_A = z1_A + c1(a) * y;
end
%% 
fourthOrderOut = smpparsvffiltbidirectional2(Bmfwd, Bmbwd, Am, double(Y)) .* divareaSum;
fourthOrderOut = fourthOrderOut(1 : halfLen);
%% plot
figure(2)
semilogx(f,20*log10(abs(Y(1 : halfLen))),':',f,20*log10(abs(Z1(1 : halfLen))))
hold on;
semilogx(f,20*log10(abs(firstOrderOut(1 : halfLen))))
semilogx(f,20*log10(abs(secondOrderOut(1 : halfLen))))
semilogx(f,20*log10(abs(fourthOrderOut)))
hold off;
grid on
axis tight;
hold off;
grid on
axis tight;
legend('Unsmoothed', 'Matrix multiplication(Gaussian window)', 'First order IIR', 'Second order IIR', 'Fourth order IIR')
title('Frequency dependent Gaussian window')
function x_oct = smoothSpectrum2(X,f,Noct)
% calculates a Gaussian function for each frequency, deriving a
% bandwidth for that frequency
x_oct = X; % initial spectrum
if Noct > 0 % don't bother if no smoothing
    for i = find(f>0,1,'first'):length(f)
        sigma = (f(i) / Noct) / pi; % standard deviation
        g = exp(-( ( (f - f(i) ) .^ 2 ) / (2 * (sigma ^ 2) ) ) ); % Gaussian
        g = g ./ sum(g);
        x_oct(i) = sum(g.*X); % calculate smoothed spectral coefficient
    end
    % remove undershoot when X is positive
    if all(X>=0)
        x_oct(x_oct<0) = 0;
    end
end
end
function Noct = estimate1NOctFromQ(Q,NFFT,db)
tau = Q * (2*pi^2) / NFFT / 0.8;
% Hyperbola : (y-max_frac_pi*pi)(y-x) = -softness, that saturates the tau to a maximum value
softness = 0.005;
max_frac_pi = 0.6;
% The hyperbola solutions are calculated
rr = roots([1,-(tau+max_frac_pi*pi),pi*tau*max_frac_pi-softness]);
% the smallest solution is the desired value
tau_h = min(rr);
tau = tau_h;
% The pole is calculated for the saturated Q
w = 1/10^(db/20);
pol = [2*w-(1+cos(tau)) -(4*w*cos(tau)-2*(1+cos(tau))) 2*w-(1+cos(tau))];
r = roots(pol);
% The pole inside the unit circle is the the desired solution
p = r(abs(r)<1);
Noct = -(3*NFFT*0.25 .* log(p))./10;
end
function y = gaussmf(x, params)
sig = params(1);
c = params(2);
y_val = @(x_val) exp((-(x_val - c)^2)/(2 * sig^2));
y = arrayfun(y_val, x);
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