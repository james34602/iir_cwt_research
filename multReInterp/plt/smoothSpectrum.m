addpath('D:\DSP_research\time-frequency-analysis\iir_cqt_research\multResByIIRInFreq')
%%
z1_A = 0.0;
z2_A = 0.0;
sigLen = 3001;
centre = 1500;
kDelta = zeros(sigLen, 1);
kDelta(centre + 1) = 1;
% oo = zeros(size(kDelta, 1), 1);
% c1 = a(2) + 2;
% c2 = (1 + a(2) + a(3)) / c1;
% d1 = (2 * b0 + b1) / c1;
% d2 = (b0 + b1 + b2) / (c1 * c2);
% for idx = 1 : size(kDelta, 1)
%     x = kDelta(idx);
%     y = x - z1_A - z2_A;
%     oo(idx) = b0 * y + d1 * z1_A + d2 * z2_A;
%     z2_A = z2_A + c2 * z1_A;
%     z1_A = z1_A + c1 * y;
% end
%%
sigma = 200;
x = 0:1:(sigLen-1);
standardGauss = gaussmf(x, [sigma, centre]);
standardGauss = standardGauss ./ sum(standardGauss);
rng(3)
sigmas = [200, 300, 200, 100, 40, 2];
[b, a, c1, c2] = gauss_precompute(sigmas);
%%
sel = 1;
dcSignal = ones(10000, 1);
% dcSignal = kDelta;
%  dcSignal = randn(10000, 1);
oo1 = zeros(size(dcSignal, 1), 1);
oo2 = zeros(size(dcSignal, 1), 1);
z0 = 0; z1 = 0;
s1 = 0; s2 = 0;
for idx = 1 : size(dcSignal, 1)
    if idx == 1000 || idx == 2000 || idx == 2700 || idx == 4500 || idx == 6000
        sel = sel + 1;
    end
    %% Coeff update
    x = dcSignal(idx);
    %% Simplified SVF for diagram
    %     y = x - z1_A - z2_A;
    %     oo1(idx) = b0(sel) * y + 2 * c2(sel) * z1_A + z2_A;
    %     z2_A = z2_A + c2(sel) * z1_A;
    %     z1_A = z1_A + c1(sel) * y;
    %% Simplified SVF for runtime
    y = x - z1_A - z2_A;
    st = c2(sel) * z1_A;
    oo1(idx) = b(sel) * y + 2 * st + z2_A;
    z2_A = z2_A + st;
    z1_A = z1_A + c1(sel) * y;
    %% Direct form
    Yi = b(sel) * x + z0;
    z0 = z1 - a(sel, 2) * Yi;
    z1 = -a(sel, 3) * Yi;
    oo2(idx) = Yi;
end
sel = 1;
figure(1)
plot(oo1);
hold on
plot(oo2);
hold off
axis tight;
legend('State variable', 'Direct form')
vYSignal = filter(b(sel), a(sel, :), kDelta);
vYSignal = filter(b(sel), a(sel, :), vYSignal(end:-1:1));
vYSignal = vYSignal(end:-1:1);
% plot(standardGauss)
% hold on;
% plot(vYSignal);
% hold off;
% axis tight;
% legend('True Gaussian', 'ARMA approximation');
%% Gaussian sigma to time constant
% fc3db = sqrt(log(2)) ./ (2 * sigmas * pi);
% alpha = cos(fc3db) - 1 + sqrt(cos(fc3db).^2 - 4*cos(fc3db) + 3);
% time = (1 - alpha) ./ alpha;
%% 1/N octave spectrum smoothing
clear
load handel.mat
Fs = 48000;
% y = zeros(2048, 1);
% y(2109) = 1;
y = y(1 : 4096);
NFFT = length(y);
if mod(NFFT,2)==0 % keep only meaningful frequencies
    halfLen = (NFFT/2)+1;
else
    halfLen = (NFFT+1)/2;
end
Y = fft(fftshift(y .* hann(NFFT)));
 Y = abs(Y);
% number of points of pre and post padding used to set initial conditions
prepad = 10;
pospad = 100;
x_fft = [Y(end-prepad+1:end); Y(1 : halfLen + pospad - 1)];
Y = Y(1:end/2+pospad);
f = ((0:halfLen-1)'./NFFT).*Fs;
thetas1 = (0:(NFFT/2+pospad)-1);
thetas1(1) = eps;
thetas1 = [fliplr(thetas1(2 : prepad + 1)), thetas1];
thetas1 = thetas1(1:NFFT/2+prepad+1);
thetas1 = [thetas1, thetas1(length(thetas1)-1:-1:1)];
thetas1 = thetas1(1 : halfLen + prepad + pospad - 1);
% frequency bins grid (linear in this case) - pre and pos padding is added
Q = 4;
% slope = pi / NFFT;
% poles = calculate_poles(slope*thetas1,Q,NFFT);
Noct = estimate1NOctFromQ(Q + (0.9/3),NFFT,3);
% thetas1 = [thetas1(end-prepad+1:end), thetas1(1 : halfLen + pospad - 1)];
%% Smoothing by time varying filter
% Noct = 8;
sigmas = (thetas1 ./ NFFT) ./ Noct / pi * NFFT;
[b, ~, c1, c2] = gauss_precompute(sigmas);
Z1 = smoothSpectrum2(Y, ((0:(NFFT/2+pospad)-1)'./NFFT).*Fs, Noct);
%% First order time varying filtering
fgt_facM = 0.3 / Noct;
poles = exp(-1 ./ (fgt_facM .* thetas1'));
tmp1 = zeros(length(x_fft), 1);
firstOrderOut = zeros(length(x_fft), 1);
% Forward filtering
z = 0;
for n=1:length(x_fft)
    o = x_fft(n) + z;
    z = poles(n) * o;
    tmp1(n) = o;
end
% Reverse filtering
z = 0;
for n=length(x_fft):-1:1
    o = tmp1(n) + z;
    z = poles(n) * o;
    firstOrderOut(n) = o;
end
firstOrderOut = firstOrderOut .* (1 - poles).^2;
%% Second order time varying filtering
tmp2 = zeros(length(x_fft), 1);
secondOrderOut = zeros(length(x_fft), 1);
% Forward filtering
z1_A = 0; z2_A = 0;
for a = 1:length(x_fft)
    x = x_fft(a);
    y = x - z1_A - z2_A;
    st = c2(a) * z1_A;
    tmp2(a) = b(a) * y + 2 * st + z2_A;
    z2_A = z2_A + st;
    z1_A = z1_A + c1(a) * y;
end
% Reverse filtering
z1_A = 0; z2_A = 0;
for a = length(x_fft):-1:1
    x = tmp2(a);
    y = x - z1_A - z2_A;
    st = c2(a) * z1_A;
    secondOrderOut(a) = b(a) * y + 2 * st + z2_A;
    z2_A = z2_A + st;
    z1_A = z1_A + c1(a) * y;
end
%% plot
figure(2)
semilogx(f,20*log10(abs(Y(1 : halfLen))),':',f,20*log10(abs(Z1(1 : halfLen))))
hold on;
semilogx(f,20*log10(abs(firstOrderOut(prepad + 1 : end - pospad + 1))))
semilogx(f,20*log10(abs(secondOrderOut(prepad + 1 : end - pospad + 1))))
hold off;
grid on
axis tight;
hold off;
grid on
axis tight;
legend('Unsmoothed', 'Matrix multiplication(Gaussian window)', 'First order IIR', 'Second order IIR')
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
function poles = calculate_poles(thetas,Q,NFFT,slope)
% function poles = calculate_poles(thetas,Q,NFFT,slope)
%
% Function that calculates the poles of a IIR-LTV-Q-FFT transform.
% For each digital frequency theta a pole is designed that complies with
% the specified Q factor, for the given number of FFT point.
%
% thetas - digital frequency grid (only half the spectrum - size: NFFT/2 + padding)
%      Q - quality factor
%   NFFT - number of points of the FFT
if nargin < 4; slope = 0.9/3; end;
N = length(thetas);
poles = zeros(N,1);
for i=1:N
    poles(i) = design_pole(thetas(i),Q+slope*thetas(i),NFFT,3);
end
end
function [p] = design_pole(theta,Q,NFFT,db)
% function p = design_pole(theta,Q,NFFT,db)
%
% Function that designs a pole for the IIR-LTV-Q-FFT transform.
% The Q factor sets the number of cycles of a sinusoid of digital frequency
% theta that must fit into the window. Thus the window is designed so that
% Q cycles fit within its db drop points (normally 3dB).
tau = Q * (2*pi^2) / NFFT / theta;
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
end