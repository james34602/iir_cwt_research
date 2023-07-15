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
% fs = 48000;
% data = [zeros(128, 1); ones(3000, 1); zeros(3000, 1)];
% Q = 16;
% fmin = 1;
% alpha = 1 + 1 / Q;
% fk = fmin * (alpha) .^ (0 : 2048);
o = emaFFTSpec_sum_(data,fs,4096,128,16,4); % Hopsize being too small causing sidelobe overwhelm the spectrum
magSpec = abs(o)';
imagesc(log10(magSpec + 0.001))
colormap(jet);
% ylim([1 50])
set(gca,'YDir','normal');
function o = emaFFTSpec_sum_(y,fs,windowSize, hopSize, Q, order)
lennon = windowSize/2+1;
freqs = 0 : fs / windowSize : fs / 2;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = 1 ./ warper * order;
firstIdx = find(time <= (windowSize / hopSize), 1, 'first');
time(1 : firstIdx) = windowSize / hopSize;
fdl = fractionalDL((max(time) - time)/2);
st = initCplxMovAvg(uint32(time));
temp = [zeros(windowSize/2,1); y; zeros(windowSize/2-1,1)];

tsAcc = 0;
o = zeros(ceil(length(y)./hopSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize - hopSize * 2) / 2;
window1(len+1:len+hopSize*2) = idealHann(hopSize * 2);
c = 1;
for a = 1:hopSize:length(temp)-windowSize
    frame = temp(a:a+windowSize-1);
    windowedFrame = frame .* window1;
    wndowedFrameShifted = circshift(windowedFrame, tsAcc);
    tsAcc = tsAcc + hopSize;
    if tsAcc >= windowSize
        tsAcc = tsAcc - windowSize;
    end
    temp2 = fft(wndowedFrameShifted);
    temp2 = temp2(1 : lennon);
    summed = procCplxMovAvg(st, temp2);
    o(c, :) = fdl.process(summed);
    c = c+1;
end
end