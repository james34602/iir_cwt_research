addpath('D:\DSP_research\time-frequency-analysis\iir_cqt_research\multResByIIRInFreq')
frameSize = 4096;
halfLen = frameSize / 2 + 1;
hop = 256;
oct = 32;
HFSamplingLimit = 0.7;
order = 2;
reqSynthesisWnd = 1;
fs = 48000;
ltv_plot(frameSize, hop, fs, oct, order, HFSamplingLimit, 0, 1, reqSynthesisWnd);
function ltv_plot(fftLen, hop, fs, oct, order, HFSamplingLimit, gaussSigmaLimiting, cornerSmoothing, reqSynthesisWnd)
zp = 1;
halfLen = fftLen / 2 + 1;
% number of points of pre and post padding used to set initial conditions
prepad = 10;
pospad = 100;
% frequency bins grid (linear in this case) - pre and pos padding is added
% poles of the IIR LTV Q FFT transform for the parameters above
f = (0:1:fftLen/2)*fs/fftLen;
% number of points of pre and post padding used to set initial conditions
thetas1 = (0:(fftLen/2+pospad)-1);
thetas1(1) = eps;
thetas1 = [fliplr(thetas1(2 : prepad + 1)), thetas1];
thetas1 = thetas1(1:fftLen/2+prepad+1);
thetas1 = [thetas1, thetas1(length(thetas1)-1:-1:1)];
thetas1 = thetas1(1 : halfLen + prepad + pospad - 1);
% thetas1 = [thetas1(end-prepad+1:end), thetas1(1 : halfLen + pospad - 1)];
%% Compute poles
if order == 2
    sigmas = (thetas1 ./ fftLen) ./ oct / pi * fftLen;
    [b, a, c1, c2] = gauss_precompute(sigmas);
else
    a = 1 - exp(-1 ./ (0.3 / oct .* thetas1.'));
    b = [];
    c1 = [];
    c2 = [];
end
%% Pole limiting
truWnd = hann(fftLen, 'periodic');
chopedWnd = truWnd(halfLen : end);
div = round(hop / (1 + HFSamplingLimit));
halfWndLen = halfLen - 1;
digw = linspace(0, pi - pi / halfWndLen, halfWndLen);
digw(halfWndLen) = pi - pi / halfWndLen;
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
s = cplxFreq(div);
if order == 2
    bDeflated = s .* s .* b(:);
    aDeflated = s .* (s + a(:, 2)) + a(:, 3);
else
    bDeflated = s .* a(:);
    aDeflated = s - (1 - a(:));
end
hHopPtCplx = bDeflated ./ aDeflated;
hHopPt = hHopPtCplx .* conj(hHopPtCplx) * chopedWnd(div);
if reqSynthesisWnd
    hHopPt = hHopPt * chopedWnd(div);
end
firstUndersampling = find((hHopPt(prepad + 1 : end - pospad + 1) * 2) <= 1, 1, 'first');
if ~isempty(prepad + firstUndersampling - 1)
    thetas1(prepad + firstUndersampling - 1 : halfLen) = thetas1(prepad + firstUndersampling - 1);
    thetas1(1 : prepad) = thetas1(prepad * 2 + 1 : -1 : prepad + 2);
    thetas1(halfLen + 1 : end) = thetas1(halfLen-1:-1:halfLen - prepad - pospad + 1);
    if cornerSmoothing
        time = zp * 0.026 * fftLen;
        alpha = 1 / (1 + time);
        [b2, a2] = butter(1, alpha);
        thetas1 = filtfilt(b2, a2, thetas1);
    end
    %% Recompute poles
    if gaussSigmaLimiting
        if order == 2
            sigmas = (thetas1 ./ fftLen) ./ oct / pi * fftLen;
            [b, a, c1, c2] = gauss_precompute(sigmas);
        else
            a = 1 - exp(-1 ./ (0.3 / oct .* thetas1.'));
        end
    end
end
if abs(thetas1 - mean(thetas1)) < eps
    flat = 1;
else
    flat = 0;
end
disp('First constant Q frequnecy that undersampling ' + string(f(firstUndersampling)) + ' Hz')
%% Compute reassignment frequency domain sample shifter
phaseShifter1 = exp(1i * (2*pi*(0:halfLen-1)/fftLen).');
phaseShifter2 = exp(-1i * (2*pi/fftLen).');
phaseShifter3 = exp(1i * (2*pi/fftLen).');
%% Obtain filterbank overall response
corrF = [];
for i = 1 : fftLen / hop - 1
    stepSize = mod(fftLen - hop * i - fftLen / 2 - 1, fftLen);
    phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
    q_fft_frameinv = ltv_1st_2ndNoSq2(phaseShifter, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen);
    %% Inverse transform
    % Virtually multiplying Hann window in time domain on frequency domain
    if reqSynthesisWnd
        q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
    end
    % Reflect spectrum and conjugate
    q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    yInvShifted = ifftshift(ifft(q_fft_frameinv)).';
    tmp(i, :) = yInvShifted;
end
n = ceil(size(tmp, 2) / (fftLen / hop));
ind = bsxfun(@plus, 1:size(tmp, 2), (0:size(tmp, 1)-1).'*n);
systemImpulse = accumarray(ind(:), tmp(:));
systemImpulse = systemImpulse(1 : fftLen - hop);
if hop ~= fftLen / 2
    truncatedSystemImpulse = [systemImpulse((fftLen - hop) - (fftLen/2) : end)];
else
    truncatedSystemImpulse = [0; systemImpulse];
end
truncatedSystemImpulse = [truncatedSystemImpulse; truncatedSystemImpulse(halfLen - 1: -1 : 1)];
% transfer = abs(fft(truncatedSystemImpulse));
% transfer = transfer(1 : length(transfer) / 2 + 1);
% transfer = interp1(linspace(0, 1, length(transfer)), transfer, linspace(0, 1, halfLen));
truncatedSystemImpulse(1) = [];
truncatedSystemImpulse = fft(truncatedSystemImpulse);
truncatedSystemImpulse = abs(truncatedSystemImpulse(1 : halfLen));
if flat
    truncatedSystemImpulse(:) = mean(truncatedSystemImpulse);
end
corrF = 1 ./ truncatedSystemImpulse;
%% Obtain unnormalized DFT filterbank frequency response and simulate result after overlap-add
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd.';
if reqSynthesisWnd
    h2 = h2 .* chopedWnd.';
end
theoreticalWindowShape = [zeros(size(thetas1, 2), 1), h2(:, (halfLen-1):-1:2), h2];
theoreticalWindowShape = theoreticalWindowShape(prepad + 1 : end - pospad + 1, :);
theoreticalWindowShape = theoreticalWindowShape .* corrF;

alt = (-1).^(0:(halfLen - 1)); % fftshift in frequency domain

dftMtx = dftmtx(fftLen);
dftMtx = dftMtx(1 : halfLen, :);
w1 = alt .* (theoreticalWindowShape .* dftMtx)';
cpxRes = fft(w1);
cpxRes = cpxRes(1 : halfLen, :);
dftFilterbank = abs(cpxRes);
figure(1)
plot(f, dftFilterbank)
axis tight
title('Frequency response of DFT filterbank')
summ = zeros(size(theoreticalWindowShape, 1), size(theoreticalWindowShape, 2));
for idx = 1 : fftLen / hop
    shifted = circshift(theoreticalWindowShape, (idx - 1) * hop, 2);
    summ = summ + shifted;
end
figure(2)
plot(summ');
axis tight
title('Frequency dependent window function after overlap-add reconstruction')
%% Obtain DFT filterbank frequency response and simulate result after overlap-add
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd.';
if reqSynthesisWnd
    h2 = h2 .* chopedWnd.';
end
theoreticalWindowShape = [zeros(size(thetas1, 2), 1), h2(:, (halfLen-1):-1:2), h2];
theoreticalWindowShape = theoreticalWindowShape(prepad + 1 : end - pospad + 1, :);
figure(3)
plot(theoreticalWindowShape')
axis tight
dftMtx = dftmtx(fftLen);
dftMtx = dftMtx(1 : halfLen, :);
w2 = alt .* (theoreticalWindowShape .* dftMtx)';
cpxRes = fft(w2);
cpxRes = cpxRes(1 : halfLen, :);
dftFilterbank = abs(cpxRes);
figure(4)
plot(f, dftFilterbank)
axis tight
title('Frequency response of DFT filterbank')
summ = zeros(size(theoreticalWindowShape, 1), size(theoreticalWindowShape, 2));
for idx = 1 : fftLen / hop
    shifted = circshift(theoreticalWindowShape, (idx - 1) * hop, 2);
    summ = summ + shifted;
end
figure(5)
plot(summ');
axis tight
title('Frequency dependent window function after overlap-add reconstruction')
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen)
%% Hann in frequency domain with input being shifted by fftshift
specHann = 2 * dftSpec + [conj(dftSpec(2)); dftSpec(1 : end - 1)] + [dftSpec(2 : end); conj(dftSpec(end - 1))];
%% Hann in frequency domain with 1 sample delayed input being shifted by fftshift
shiftedspecHann = phaseShifter1 .* (2 * dftSpec + phaseShifter2 .* [conj(dftSpec(2)); dftSpec(1 : end - 1)] + phaseShifter3 .* [dftSpec(2 : end); conj(dftSpec(end - 1))]);
% Hann windowed
x_fft1 = [conj(specHann(prepad + 1 : -1 : 2, :)); specHann; conj(specHann(halfLen - 1 : -1 : (halfLen - pospad + 1), :))] / 4;
x_fft2 = [conj(shiftedspecHann(prepad + 1 : -1 : 2, :)); shiftedspecHann; conj(shiftedspecHann(halfLen - 1 : -1 : (halfLen - pospad + 1), :))] / 4;
% Rectangular windowed
% x_fft1 = [conj(dftSpec(prepad + 1 : -1 : 2, :)); dftSpec; conj(dftSpec(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
% smpShifted = dftSpec .* phaseShifter1;
% x_fft2 = [conj(smpShifted(prepad + 1 : -1 : 2, :)); smpShifted; conj(smpShifted(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
%% Gaussian windowing
x_fft = [x_fft1, x_fft2];
tmp = zeros(size(x_fft, 1), 2, 'like', x_fft);
q_fft_frame = ltv(x_fft, tmp, b, a, c1, c2);
% Remove periodic padding
spec = q_fft_frame(prepad + 1 : end - pospad + 1, 1);
if ~isempty(corrF)
    spec = spec .* corrF;
end
spec(1, :) = real(spec(1, :));
spec(end, :) = real(spec(end, :));
end