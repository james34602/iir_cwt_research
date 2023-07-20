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
fftLen = 4096;
hopSize = 128;
% t = 0:1/fs:1-1/fs;
% dt = seconds(t(2)-t(1));
% rng(1)
% x1 = chirp(t,100,t(end),20800,'quadratic');
% x2 = exp(2j*pi*3800*cos(2*pi*2*t));
% data = real((x1 + x2 + randn(size(t))/10))';
o = emaFFTSpec_m4(data, fs, fftLen, hopSize, 32, 4); % Hopsize being too small causing sidelobe overwhelm the spectrum
% o = traditional(data, fs, fftLen, hopSize);
% o = o(fftLen / hopSize + 1 : end, :);
magSpec = abs(o)';
imagesc(log10(magSpec + 0.001))
colormap(jet);
set(gca,'YDir','normal');
function o = emaFFTSpec_m4(y, fs, windowSize, hopSize, Q, order)
lennon = windowSize/2+1;
freqs = 0 : fs / windowSize : fs / 2;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = (1 ./ warper) / 2; % This is window size, impulse in the middle
% time = time - time(end);
firstIdx = find(time <= (windowSize / hopSize) / 2, 1, 'first');
time(1 : firstIdx) = (windowSize / hopSize) / 2;
% time = time / order;
if any(isinf(time))
    error('Inf detected');
end
fdl = fractionalDL(max(time) - time);
fgt_fac = 1 ./ (1 + time);
fc = 2 * (tan(log(1 - fgt_fac) / -2)) / pi;
coeffB = zeros(lennon, 1);
coeffA = zeros(lennon, order + 1);
for idx = 1 : lennon
    [coeffB(idx, :), coeffA(idx, :)] = bessel(order, time(idx));
end
coeffA = coeffA.';
st = initNPole1Zero(uint32(lennon), uint32(order));
tsAcc = 0;
o = zeros(ceil(length(y)./hopSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize - hopSize * 2) / 2;
window1(len + 1 : len + hopSize * 2) = hann(hopSize * 2, 'periodic');
temp = [zeros(windowSize / 2, 1); y; zeros(windowSize / 2 - 1, 1)];
c = 1;
% states = zeros(lennon, order);
for a = 1:hopSize:length(temp)-windowSize
    frame = circshift(temp(a:a+windowSize-1) .* window1, tsAcc);
    tsAcc = tsAcc + hopSize;
    if tsAcc >= windowSize
        tsAcc = tsAcc - windowSize;
    end
    temp2 = fft(frame);
    temp2 = temp2(1 : lennon);
%     for idx = 1 : lennon
%         [o(c, idx), states(idx, :)] = filter(coeffB(idx, :), coeffA(idx, :), temp2(idx), states(idx, :));
%     end
    filterOut = procNPole1Zero(st, coeffB, coeffA, temp2);
    o(c, :) = fdl.process(filterOut);
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