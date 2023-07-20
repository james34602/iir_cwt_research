% One pole IIR, coefficient to time in sample
% time = (1 - alpha) / alpha
%
% One pole IIR coeff, time in sample to coefficient
% alpha = 1 / (1 + time)
%
% One pole IIR coeff, frequency(-3 dB point) to coefficient
% alpha = 1 - exp(-2 * atan(pi * fcNorm / 2))
%
% One pole IIR coeff, coefficient to frequency(-3 dB point)
% fcNorm = 2 * (tan(log(1 - alpha) / -2)) / pi
% data = [zeros(128, 1); 1; zeros(8192, 1)];
fftLen = 512;
hopSize = 16;
%%
load signal1
rng(1);
data = [rand(64, 1); s; zeros(256, 1)];
fs = 48000;
[X_complex, ez, timeComp] = emaFFTSpec_m3(data, fs, fftLen, hopSize, 32, 5); % Hopsize being too small causing sidelobe overwhelm the spectrum
reassign_frequency = 1;
tfrMag = 0;
frameSize = fftLen;
halfLen = frameSize / 2 + 1;
X_complex = X_complex.';
X_mag = abs(X_complex);
delayed = circshift(data, 1);
delayed(1) = 0;
X_prev_time = emaFFTSpec_m3(delayed, fs, fftLen, hopSize, 32, 5); % Hopsize being too small causing sidelobe overwhelm the spectrum
fdl = fractionalDL(timeComp);
X_prev_time = X_prev_time.';
% cross-spectra - ie. spectra of cross-correlation between the
% respective time-domain signals
X_cross_time = cross_spectrum(X_complex, X_prev_time);
% instantaneous frequency estimates
% normalized frequencies in range [0.0, 1.0] - from DC to sample rate
X_inst_freqs = arg(X_cross_time);

X_mag = X_mag(1 : halfLen, :);
% X_mag(1, :) = X_mag(1, :) * 0.5;
X_complex = X_complex(1 : halfLen, :);
X_inst_freqs = X_inst_freqs(1 : halfLen, :);
% Reassigned spectrogram requantized both in frequency and time.
% Note it is quantized into non-overlapping output time frames which may be
% of a different size than input time frames.
% group delays are in range [-0.5, 0.5] - relative coordinates within the
% frame where 0.0 is the frame center
if reassign_frequency
    X_y = X_inst_freqs;
else
    X_y = repmat(linspace(0, 0.5, halfLen), [size(X_inst_freqs, 2), 1]).';
end
if tfrMag
    Snew = X_mag;
else
    Snew = X_complex;
end
RS = zeros(size(X_y, 1), 1);
RS2 = zeros(size(X_y));
range2 = floor(mod(X_y * size(X_y, 1) * 2, frameSize)) + 1;
range2(range2 < 1) = 1;
range2(range2 > size(X_y, 1)) = size(X_y, 1);
for k = 1 : size(X_y, 2)
    mag = Snew(:, k) .* ez(k, :).';
    RS(:) = 0;
    for j = 1 : halfLen
        RS(range2(j, k)) = RS(range2(j, k)) + mag(j);
    end
    RS2(:, k) = fdl.process(RS);
end
RS2 = RS2(:, 1 : end);
if tfrMag
    imagesc((abs(RS2) + 0.001).^alpha)
else
    imagesc((abs(RS2) + 0.001).^alpha)
end
colormap(jet);
set(gca,'YDir','normal');
function [o, sht, timeComp] = emaFFTSpec_m3(y,fs,windowSize,hopSize,Q,order)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = (1 ./ warper); % This is window size, impulse in the middle
firstIdx = find(time <= (windowSize / hopSize), 1, 'first');
time(1 : firstIdx) = (windowSize / hopSize);
if isempty(firstIdx)
    time(:) = (windowSize / hopSize);
end
time = time / 2;
timeComp = (max(time) - time);
time = time / order;
if any(isinf(time))
    error('Inf detected');
end
fgt_fac = 1 ./ (1 + time);
lennon = windowSize/2+1;
tsAcc = 0;
o = zeros(ceil(length(y)./hopSize),lennon);
sht = zeros(ceil(length(y)./hopSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize - hopSize * 2) / 2;
window1(len + 1 : len + hopSize * 2) = hann(hopSize * 2, 'periodic');
temp = [zeros(windowSize / 2, 1); y; zeros(windowSize / 2 - 1, 1)];
c = 1;
ema = zeros(lennon,order);
for a = 1:hopSize:length(temp)-windowSize
    windowedFrame = temp(a:a+windowSize-1) .* window1;
    frame = circshift(windowedFrame, tsAcc);
    oldShift = tsAcc;
    inds = 0 : (windowSize - 1);
    m = oldShift;
    ez = exp(-1i*2*pi*m.*inds./windowSize)';
    tsAcc = tsAcc + hopSize;
    if tsAcc >= windowSize
        tsAcc = tsAcc - windowSize;
    end
    temp2 = fft(frame);
    temp2 = temp2(1 : lennon);
    ema(:, 1) = ema(:, 1) + fgt_fac .* (temp2 - ema(:, 1));
    for b = 2 : order
        ema(:, b) = ema(:, b) + fgt_fac .* (ema(:, b - 1) - ema(:, b));
    end
    sht(c, :) = ez(1 : lennon);
    o(c, :) = ema(:, order);
%     yf = o(c, :);
%     yf(lennon+1:windowSize) = conj(yf(lennon-1:-1:2));
%     aa = real(ifft(yf));
%     inverseShift = circshift(aa, -oldShift);
%     figure(1)
%     plot(windowedFrame)
%     hold on;
%     plot(inverseShift);
%     hold off;
%     axis tight;
%     o(c, :) = fdl.process(ema(:, order));
    c = c+1;
end
end
function o = traditional(y, fs, windowSize, hopSize)
lennon = windowSize/2+1;
temp = [zeros(windowSize,1); y; zeros(windowSize,1)];
o = zeros(ceil(length(y)./hopSize),lennon);
window1 = hann(windowSize, 'periodic');
c = 1;
for a = 1:hopSize:length(temp)-windowSize
    frame = temp(a:a+windowSize-1) .* window1;
    temp2 = fft(frame);
    temp2 = temp2(1 : lennon);
    o(c, :) = temp2;
    c = c+1;
end
end
function y = cross_spectrum(spectrumA, spectrumB)
%     Returns a cross-spectrum, ie. spectrum of cross-correlation of two signals.
%     This result does not depend on the order of the arguments.
%     Since we already have the spectra of signals A and B and and want the
%     spectrum of their cross-correlation, we can replace convolution in time
%     domain with multiplication in frequency domain.
y = spectrumA .* conj(spectrumB);
end
function y = arg(values)
%     Argument (angle) of complex numbers wrapped and scaled to [0.0, 1.0].
%
%     input: an array of complex numbers
%     output: an array of real numbers of the same shape
%
%     np.angle() returns values in range [-np.pi, np.pi].
y = mod(angle(values) / (2 * pi), 1.0);
end