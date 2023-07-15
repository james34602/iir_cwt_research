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
o = emaFFTSpec_m3(data, fs, fftLen, hopSize, 32, 5); % Hopsize being too small causing sidelobe overwhelm the spectrum
o2 = traditional(data, fs, fftLen, hopSize);
figure(1)
magSpec = abs(o)';
imagesc(log10(magSpec + 0.001))
colormap(jet);
set(gca,'YDir','normal');
ylim([1, 50])
figure(2)
magSpec = abs(o2)';
imagesc(log10(magSpec + 0.001))
colormap(jet);
set(gca,'YDir','normal');
ylim([1, 50])
function o = emaFFTSpec_m3(y,fs,windowSize,hopSize,Q,order)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = (1 ./ warper); % This is window size, impulse in the middle
firstIdx = find(time <= (windowSize / hopSize), 1, 'first');
time(1 : firstIdx) = (windowSize / hopSize);
time = time / 2;
fdl = fractionalDL((max(time) - time));
time = time / order;
if any(isinf(time))
    error('Inf detected');
end
fgt_fac = 1 ./ (1 + time);
lennon = windowSize/2+1;
tsAcc = 0;
o = zeros(ceil(length(y)./hopSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize - hopSize * 2) / 2;
window1(len + 1 : len + hopSize * 2) = hann(hopSize * 2, 'periodic');
temp = [zeros(windowSize / 2, 1); y; zeros(windowSize / 2 - 1, 1)];
c = 1;
% time = 1;
% order = 5;
% kDelta = [1; zeros(1023, 1)];
% yy = zeros(order, 1024);
% legStr = [];
% for idx = 1 : order
%     a = 1 ./ (1 + time / idx);
%     lst = [];
%     for ii = 1 : idx
%         lst = [lst, dfilt.df1(a, [1, -(1 - a)])];
%     end
%     yy(idx, :) = filter(dfilt.cascade(lst), kDelta);
%     legStr = [legStr, "Order" + string(idx)];
% end
% plot(yy');
% axis tight;
% legend(legStr)
ema = zeros(lennon,order);
for a = 1:hopSize:length(temp)-windowSize
    windowedFrame = temp(a:a+windowSize-1) .* window1;
    frame = circshift(windowedFrame, tsAcc);
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
    o(c, :) = fdl.process(ema(:, order));
    c = c+1;
end
end
function o = traditional(y, fs, windowSize, hopSize)
lennon = windowSize/2+1;
temp = [zeros(windowSize,1); y; zeros(windowSize,1)];
o = zeros(ceil(length(y)./hopSize),lennon);
window1 = idealHann(windowSize);
c = 1;
for a = 1:hopSize:length(temp)-windowSize
    frame = temp(a:a+windowSize-1) .* window1;
    temp2 = fft(frame);
    temp2 = temp2(1 : lennon);
    o(c, :) = temp2;
    c = c+1;
end
end