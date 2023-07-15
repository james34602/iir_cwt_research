frameSize = 4096;
fs = 48000;
oct = 28;
wnd = generateGaussianWindowOnFreqDomain(frameSize,fs,oct);
function wnd = generateGaussianWindowOnFreqDomain(fftLen,fs,oct)
halfLen = fftLen / 2 + 1;
f = (0:1:fftLen/2)*fs/fftLen;
%% Ideal Gaussian
freqCentreDiff = (repmat(f, halfLen, 1)' - f ) .^ 2;
flat = fft(fftshift(hann(fftLen, 'periodic')));
flat = real(flat(1 : halfLen));
wnd = zeros(halfLen, halfLen);
sigmas = f ./ oct / pi;
sigmas(sigmas == 0) = eps;
sigmas = -1 ./ (2 * (sigmas .^ 2));
for j = 1:halfLen
    g2 = exp(freqCentreDiff .* sigmas(j)); % Gaussian
    transformationMatrix = g2 ./ sum(g2, 2);
    res = transformationMatrix * flat;
    wnd(j, :) = res;
end
wnd(:, halfLen+1:fftLen) = conj(wnd(:, halfLen-1:-1:2));
wnd = ifftshift(ifft(wnd.'), 1);
wnd = wnd - min(wnd);
wnd = wnd ./ max(wnd);
end