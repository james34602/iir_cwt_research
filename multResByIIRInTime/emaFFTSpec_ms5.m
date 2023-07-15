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
fs = 48000;
data = [zeros(512, 1); 1; zeros(16384, 1)];
fftLen = 4096;
hopSize = 512;
dbg = 0;
[o, noMultRes, multRes] = emaFFTSpec_m3(data, fs, fftLen, hopSize, 32, 5, 0, 0); % Hopsize being too small causing sidelobe overwhelm the spectrum
% plot(data);
% hold on
% plot(noMultRes(hopSize + 1 : end))
plot(multRes(hopSize + 1 : end))
% hold off
axis tight
ylim([-0.005 0.01])
% magSpec = abs(o)';
% imagesc(log10(magSpec + 0.001))
% colormap(jet);
% set(gca,'YDir','normal');
function [o, noMultRes, multRes] = emaFFTSpec_m3(y,fs,windowSize,hopSize,Q,order,debug1,debug2)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = (1 ./ warper); % This is window size, impulse in the middle
firstIdx = find(time <= (windowSize / hopSize), 1, 'first');
time(1 : firstIdx) = (windowSize / hopSize);
time = time / 2.1;
mm = max(time) - time;
fdl = fractionalDL(mm);
time = time / order;
if any(isinf(time))
    error('Inf detected');
end
fgt_fac = 1 ./ (1 + time);
% fgt_fac(:) = fgt_fac(1);
% fgt_fac(:) = 1;
a = fgt_fac(1);
% lst = [];
% for ii = 1 : order
%     lst = [lst, dfilt.df1(a, [1, -(1 - a)])];
% end
% fvtool(dfilt.cascade(lst));
lennon = windowSize / 2 + 1;
tsAcc = 0;
o = zeros(ceil(length(y) / hopSize), lennon);
window1 = zeros(windowSize, 1);
len = (windowSize - hopSize * 2) / 2;
window1(len + 1 : len + hopSize * 2) = hann(hopSize * 2, 'periodic');
% window1(len + 1 : len + hopSize * 2) = 1;
temp = [zeros(windowSize / 2, 1); y; zeros(windowSize / 2 - 1, 1)];
window2 = zeros(windowSize,1);
window2(len + 1 : len + hopSize * 2) = 1;
c = 1;
ema = zeros(lennon,order);
noMultRes = zeros(ceil(length(y) / hopSize) * hopSize, 1);
multRes = zeros(ceil(length(y) / hopSize) * hopSize, 1);
accFrame1 = zeros(windowSize, 1);
accFrame2 = zeros(windowSize, 1);
for a = 1:hopSize:length(temp)-windowSize
    originalFrame = temp(a:a+windowSize-1) .* window2;
    windowedFrame = temp(a:a+windowSize-1) .* window1;
    frame = circshift(windowedFrame, tsAcc);
    oldShift = tsAcc;
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
    dline = fdl.process(ema(:, end));
    o(c, :) = ema(:, end);
    %% With multiresolution
    yf = o(c, :).';
    yf(lennon+1:windowSize) = conj(yf(lennon-1:-1:2));
    aa = real(ifft(yf));
    inverseShift = circshift(aa, -oldShift - len - hopSize * 2);
    accFrame1 = accFrame1 + inverseShift;
    Ss = hopSize;
    myOut = accFrame1(1 : Ss);
    accFrame1 = [accFrame1(Ss + 1 : end); zeros(Ss, 1)];
    if debug1
        figure(1)
        plot(inverseShift)
        hold on;
        plot(accFrame1);
        hold off;
        axis tight;
        figure(2)
        magSpec = abs(o)';
        imagesc(log10(magSpec + 0.001))
        colormap(jet);
        set(gca,'YDir','normal');
    end
    multRes(a:a+hopSize-1) = myOut;
    %% No multiresolution
    yf = temp2;
    yf(lennon+1:windowSize) = conj(yf(lennon-1:-1:2));
    aa = ifft(yf);
    inverseShift = circshift(aa, -oldShift);
    accFrame2 = accFrame2 + inverseShift;
    myOut = accFrame2(1 : Ss);
    accFrame2 = [accFrame2(Ss + 1 : end); zeros(Ss, 1)];
    noMultRes(a:a+hopSize-1) = myOut;
    if debug2
%         plot(originalFrame(len + 1 : len + hopSize));
%         hold on
        plot(accFrame2);
%         hold off;
    end
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