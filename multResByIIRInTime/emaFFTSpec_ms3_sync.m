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
hopSize = 1;
%%
load signal1
rng(1);
s = [rand(64, 1); s; zeros(300, 1)];
%%
data = s;
[o, o2] = emaFFTSpec_m3(data, fs, fftLen, hopSize, 12, 5); % Hopsize being too small causing sidelobe overwhelm the spectrum
magSpec = abs(o2)';
alpha = 0.15;             %% TFR compression factor (modify display contrast)
imagesc((magSpec + 0.001).^alpha)
colormap(jet);
set(gca,'YDir','normal');
function [o, o2] = emaFFTSpec_m3(y,fs,windowSize,hopSize,Q,order)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
%% https://www.desmos.com/calculator/017zevkgfk
warper = (freqs' ./ fs) ./ Q .* hopSize;
% time = exp(-warper) ./ (1 - exp(-warper)) * order;
time = (1 ./ warper); % This is window size, impulse in the middle
firstIdx = find(time <= (windowSize / hopSize), 1, 'first');
time(1 : firstIdx) = (windowSize / hopSize);
time = time / 2;
% time(:) = 16;
nonSyncTime = time + 1; % One more DFT frame
fdl1 = fractionalDL((max(time) - time));
fdl2 = fractionalDL((max(time) - time));
time = time / order;
nonSyncTime = nonSyncTime / order;
if any(isinf(time))
    error('Inf detected');
end
fgt_fac = 1 ./ (1 + time);
fgt_facNonSync = 1 ./ (1 + nonSyncTime);
lennon = windowSize/2+1;
tsAcc = 0;
o = zeros(ceil(length(y)./hopSize),lennon);
o2 = zeros(ceil(length(y)./hopSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize - hopSize * 2) / 2;
window1(len + 1 : len + hopSize * 2) = hann(hopSize * 2, 'periodic');
temp = [zeros(windowSize / 2, 1); y; zeros(windowSize / 2 - 1, 1)];
c = 1;
m = 1 : lennon;
lambda = (m'-1) ./ windowSize;
pTs    = -1.0./(nonSyncTime * 2 * hopSize) + 1i*2*pi*lambda;
emaRegular = zeros(lennon, order);
emaSync = zeros(lennon, order);
synced = zeros(lennon, 1);
inds = 0 : (windowSize - 1);
m = floor((nonSyncTime * 2 * hopSize) / 2)';
ez = exp(-1i*2*pi.*m.*inds(1 : lennon)./(nonSyncTime' * 2 * hopSize))';
fout=((0:(windowSize - 1)) ./ windowSize * 2 * pi)';
fmax = fout(end);
fout = fout(1 : lennon);
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
    emaRegular(:, 1) = emaRegular(:, 1) + fgt_facNonSync .* (temp2 - emaRegular(:, 1));
    emaSync(:, 1) = emaSync(:, 1) + fgt_fac .* (temp2 - emaSync(:, 1));
    for b = 2 : order
        emaRegular(:, b) = emaRegular(:, b) + fgt_facNonSync .* (emaRegular(:, b - 1) - emaRegular(:, b));
        emaSync(:, b) = emaSync(:, b) + fgt_fac .* (emaSync(:, b - 1) - emaSync(:, b));
    end
%     m = floor(oldShift / 2);
%     ez = exp(-1i*2*pi*m*inds/windowSize)';

    synced(:) = 0;
    for n = 1 : lennon
        if abs(emaRegular(n, order)) ~= 0
            Ts_yDg = emaSync(n, order) + pTs(n) * emaRegular(n, order);
            mrhat = round(windowSize * imag(Ts_yDg / emaRegular(n, order)) / (2.0 * pi));
            if (mrhat > windowSize)
                mrhat = mrhat - windowSize;
            end
            if (mrhat < 1)
                synced(1) = synced(1) + nonlin(emaRegular(n, order));
            elseif (mrhat >= lennon)
                synced(lennon) = synced(lennon) + nonlin(emaRegular(n, order));
            else
                synced(mrhat + 1) = synced(mrhat + 1) + nonlin(emaRegular(n, order));
            end
        else
            synced(n) = nonlin(emaRegular(n, order)); % Overwrite the result
        end
    end
    o(c, :) = fdl1.process(emaRegular(:, order));
    o2(c, :) = fdl2.process(synced);
    c = c+1;
end
end
function y = nonlin(x)
y = (x);
end