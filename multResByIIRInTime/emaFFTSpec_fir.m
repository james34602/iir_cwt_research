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
fftLen = 4096;
halfLen = fftLen / 2 + 1;
[o, time] = emaFFTSpec_fir2(data,fs,fftLen, 128, 32); % Hopsize being too small causing sidelobe overwhelm the spectrum
% time = time - time(end);
delay = zeros(size(time));
oCell = num2cell(o, 1);
wndSet = cell(1, halfLen);
for idx = 1 : halfLen
    wndSet{idx} = hann(time(idx) * 2, 'periodic');
    if isempty(wndSet{idx})
        wndSet{idx} = 1;
    end
    delay(idx) = length(wndSet{idx}) / 2;
end
compensation = floor(max(delay) - delay);
for idx = 1 : halfLen
    wndSet{idx} = [zeros(compensation(idx), 1); wndSet{idx}];
end
func = @(fir, x) filter(fir, 1, x);
o2 = cellfun(func, wndSet, oCell, 'un', 0);
o2 = cell2mat(o2);
figure(1)
magSpec = abs(o)';
imagesc(log10(magSpec + 0.001))
colormap(jet);
set(gca,'YDir','normal');
figure(2)
magSpec = abs(o2)';
imagesc(log10(magSpec + 0.001))
colormap(jet);
set(gca,'YDir','normal');
function [o, time] = emaFFTSpec_fir2(y, fs, windowSize, hopSize, Q)
lennon = windowSize/2+1;
freqs = 0 : fs / windowSize : fs / 2;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = 1 ./ warper;
firstIdx = find(time <= (windowSize / hopSize), 1, 'first');
time(1 : firstIdx) = windowSize / hopSize;
temp = [zeros(windowSize/2,1); y; zeros(windowSize/2-1,1)];
tsAcc = 0;
o = zeros(ceil(length(y)./hopSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize - hopSize * 2) / 2;
window1(len+1:len+hopSize*2) = hann(hopSize * 2, 'periodic');
c = 1;
for a = 1:hopSize:length(temp)-windowSize
    frame = circshift(temp(a:a+windowSize-1) .* window1, tsAcc);
    tsAcc = tsAcc + hopSize;
    if tsAcc >= windowSize
        tsAcc = tsAcc - windowSize;
    end
    temp2 = fft(frame);
    o(c, :) = temp2(1 : lennon);
    c = c+1;
end
end