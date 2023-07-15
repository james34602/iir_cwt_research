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
% data = [zeros(512, 1); 1; zeros(16384, 1)];
fftLen = 4096;
hopSize = 512;
data2 = [data; zeros(fftLen - hopSize, 1)];
dbg = 0;
[o, noMultRes, multRes] = emaFFTSpec_m3(data2, fs, fftLen, hopSize, 32, 0); % Hopsize being too small causing sidelobe overwhelm the spectrum
plot(data2);
hold on
plot(noMultRes(hopSize + 1 : end))
plot(multRes(hopSize + 1 : end))
hold off
axis tight
 ylim([-0.0004 0.001])
%  magSpec = abs(o)';
%  imagesc(log10(magSpec + 0.001))
%  colormap(jet);
%  set(gca,'YDir','normal');
function [o, noMultRes, multRes] = emaFFTSpec_m3(y,fs,windowSize,hopSize,Q,debug2)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = (1 ./ warper); % This is window size, impulse in the middle
firstIdx = find(time <= (windowSize / hopSize) - 1, 1, 'first');
time(1 : firstIdx) = (windowSize / hopSize) - 1;
time = uint32(time);
% time(:) = 8;
if any(isinf(time))
    error('Inf detected');
end
halfLen = windowSize / 2 + 1;
st = initCplxMovAvg(time);
% rec = [];
% rec = [rec, procCplxMovAvg(st, ones(halfLen, 1))];
% for ii = 1 : 32
%     rec = [rec, procCplxMovAvg(st, zeros(halfLen, 1))];
% end
% rec = real(rec);
lennon = windowSize / 2 + 1;
tsAcc = 0;
o = zeros(ceil(length(y) / hopSize), lennon);
window1 = zeros(windowSize, 1);
len = (windowSize - hopSize * 2) / 2;
window1(len + 1 : len + hopSize * 2) = idealHann(hopSize * 2);
% window1(len + 1 : len + hopSize * 2) = 1;
temp = [zeros(windowSize / 2, 1); y; zeros(windowSize / 2 - 1, 1)];
c = 1;
noMultRes = zeros(ceil(length(y) / hopSize) * hopSize, 1);
multRes = zeros(ceil(length(y) / hopSize) * hopSize, 1);
accFrame1 = zeros(windowSize, 1);
accFrame2 = zeros(windowSize, 1);
for a = 1:hopSize:length(temp)-windowSize
    windowedFrame = temp(a:a+windowSize-1) .* window1;
    frame = circshift(windowedFrame, tsAcc);
    oldShift = tsAcc;
    tsAcc = tsAcc + hopSize;
    if tsAcc >= windowSize
        tsAcc = tsAcc - windowSize;
    end
    temp2 = fft(frame);
    temp2 = temp2(1 : lennon);
    summed = procCplxMovAvg(st, temp2);
    o(c, :) = summed;
    %% With multiresolution
    yf = o(c, :).';
    yf(lennon+1:windowSize) = conj(yf(lennon-1:-1:2));
    aa = real(ifft(yf));
%     window2 = circshift(window1, -oldShift - len - hopSize * 1);
    inverseShift1 = circshift(aa, -oldShift - len - hopSize * 2) .* 1;
    accFrame1 = accFrame1 + inverseShift1;
    Ss = hopSize;
    myOut1 = accFrame1(1 : Ss);
    accFrame1 = [accFrame1(Ss + 1 : end); zeros(Ss, 1)];
    multRes(a:a+hopSize-1) = myOut1;
    %% No multiresolution
    yf = temp2;
%     yf(30 : end) = 0;
    yf(lennon+1:windowSize) = conj(yf(lennon-1:-1:2));
    bb = ifft(yf);
    inverseShift2 = circshift(bb, -oldShift);
    accFrame2 = accFrame2 + inverseShift2;
    myOut2 = accFrame2(1 : Ss);
    accFrame2 = [accFrame2(Ss + 1 : end); zeros(Ss, 1)];
    noMultRes(a:a+hopSize-1) = myOut2;
    if debug2
        figure(1)
        plot(aa)
        hold on;
        plot(inverseShift1);
        hold off;
        axis tight
        figure(2)
        plot(bb)
        hold on
        plot(inverseShift2);
        hold off;
        axis tight;
        disp('')
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