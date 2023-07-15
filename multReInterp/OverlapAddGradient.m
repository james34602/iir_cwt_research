rng(1)
fs = 48000;
frameSize = 1024;
halfLen = frameSize / 2 + 1;
hop = 256;
oct = 32;
HFSamplingLimit = 0.7;
order = 2;
reqSynthesisWnd = 1;
[coeff, f_q] = ltv_precomute(frameSize, hop, fs, oct, order, HFSamplingLimit, 1, 1, reqSynthesisWnd);
%% Generate chirp
D = frameSize / fs;            % duration in seconds
t = (0:frameSize-1) / fs;         % discrete-time axis (sec)
f1 = 1 * fs / frameSize;
f2 = halfLen * fs / frameSize;
bandlimitChirp = chirp2(t, D, f1, f2)';
[spec1, t_q] = ltv_spectrogram(bandlimitChirp, coeff);
S = spec1(:, :, 1);
nframes = size(S, 2);
halfLen = frameSize / 2 + 1;
wnd = hann(frameSize, 'periodic');
y2 = zeros((nframes + frameSize / hop - 1) * hop, 1);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S(:, i);
    q_fft_frameinvBk = q_fft_frameinv;
    q_fft_frameinv(halfLen+1:frameSize,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    correctedTime = ifft(q_fft_frameinv);
    correctedTime2 = correctedTime .* coeff.correctionWnd;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv = q_fft_frameinvBk .* coeff.wndCorrectionWeighting + getbackCorrectedToSpectrum2(1 : halfLen) .* (1 - coeff.wndCorrectionWeighting);
    % Reflect spectrum and conjugate
    q_fft_frameinv(halfLen+1:frameSize) = conj(q_fft_frameinv(halfLen-1:-1:2));
    yInvShifted = ifftshift(ifft(q_fft_frameinv)) .* wnd;
    % Overlap add
    y2(1+(i-1)*hop : frameSize+(i-1)*hop) = y2(1+(i-1)*hop : frameSize+(i-1)*hop) + yInvShifted;
end
y = y2(frameSize - hop + 1 : frameSize - hop + frameSize);
SNR = 10*log10(sum(abs(bandlimitChirp).^2)/sum(abs(bandlimitChirp-y).^2))
err = abs(fft(y - bandlimitChirp));
err = err(1 : halfLen);
[bSm, aSm] = butter(1, 0.09);
sm = filtfilt(bSm, aSm, err);
[bSm2, aSm2] = butter(1, 0.005);
sm2 = filtfilt(bSm2, aSm2, err);
sm2 = sm2 ./ max(sm2);
[~, idx] = sort(sm);
% Take the frequency from max power
transitionBin = idx(int32(halfLen));
correctionWndTransition = zeros(1, size(coeff.theoreticalWindowShape, 2));
for idx = 1 : frameSize / hop
    shifted1 = circshift(coeff.theoreticalWindowShape(transitionBin, :), (idx - 1) * hop, 2);
    correctionWndTransition = correctionWndTransition + shifted1;
end
correctionWndTransition = 1 ./ correctionWndTransition(:);
ny = frameSize - hop + frameSize + frameSize / 2;
nframes = fix(ceil(ny/hop) - frameSize / hop / 2);
halfLen = frameSize / 2 + 1;
wnd = hann(frameSize, 'periodic');
shiftedChirp1 = circshift(bandlimitChirp, hop / 2);
shiftedChirp2 = circshift(bandlimitChirp, hop / 4);
[spec1, t_q] = ltv_spectrogram(bandlimitChirp, coeff);
[spec2, t_q] = ltv_spectrogram(shiftedChirp1, coeff);
[spec3, t_q] = ltv_spectrogram(shiftedChirp2, coeff);
S1 = spec1(:, :, 1);
S1(halfLen+1:frameSize,:) = conj(S1(halfLen-1:-1:2,:));
STime1 = ifft(S1);
S2 = spec2(:, :, 1);
S2(halfLen+1:frameSize,:) = conj(S2(halfLen-1:-1:2,:));
STime2 = ifft(S2);
S3 = spec3(:, :, 1);
S3(halfLen+1:frameSize,:) = conj(S3(halfLen-1:-1:2,:));
STime3 = ifft(S3);
normalizationVector = ones(frameSize, 1) / frameSize;
%% New optimization loop 3
correctionWndHF = coeff.correctionWnd;
wndCorrectionWeightingLF = coeff.wndCorrectionWeighting;
wndCorrectionWeightingHF = 1 - wndCorrectionWeightingLF;
weights1 = [wndCorrectionWeightingLF; wndCorrectionWeightingHF];
func1 = @(x) singleWndSplitLFHF(x, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp);
[f, g] = func1(weights1);

func1_ = @(x) singleWndSplitLFHFDualCost(x, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp, S2, STime2, shiftedChirp1, S3, STime3, shiftedChirp2);
[f, g] = func1_(weights1);

weights2 = [wndCorrectionWeightingLF; wndCorrectionWeightingHF; correctionWndHF];
func2 = @(x) singleWndOptSplitLFHF(x, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp);
[f, g] = func2(weights2);
weights3 = [wndCorrectionWeightingLF; abs(randn(halfLen, 1) * 1e-8); correctionWndTransition];
func3 = @(x) dualWndSplitLFMF(x, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp);
[f, g] = func3(weights3);
weights4 = [wndCorrectionWeightingLF; abs(randn(halfLen, 1) * 1e-8); wndCorrectionWeightingHF; correctionWndTransition];
func4 = @(x) dualWndSplitLFMFHF(x, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp);
[f, g] = func4(weights4);
weightsOpt = fminunc(func1_, weights1, optimoptions('fminunc','OptimalityTolerance',1e-14,'StepTolerance',1e-14,'MaxIterations',2000,'SpecifyObjectiveGradient',true, 'Display','iter'));
function x = chirp2(t,t1,f0,f1)
beta = (f1-f0)./t1;
x = cos(2*pi * ( 0.5* beta .* (t .* t) + f0 * t));
end
function [f, g] = singleWndSplitLFHF(weights, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S1(:, i);
    correctedTime2 = STime1(:, i) .* correctionWndHF;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse = sum(meaSqr .* normalizationVector);
f = mse.value;
g = real(mse.backprop(1));
end
function [f, g] = singleWndSplitLFHFDualCost(weights, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp, S2, STime2, shiftedChirp1, S3, STime3, shiftedChirp2)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S1(:, i);
    correctedTime2 = STime1(:, i) .* correctionWndHF;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse1 = sum(meaSqr .* normalizationVector);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S2(:, i);
    correctedTime2 = STime2(:, i) .* correctionWndHF;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - shiftedChirp1) .^ 2;
mse2 = sum(meaSqr .* normalizationVector);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S3(:, i);
    correctedTime2 = STime3(:, i) .* correctionWndHF;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - shiftedChirp2) .^ 2;
mse3 = sum(meaSqr .* normalizationVector);
mse = (mse1 + mse2 + mse3) / 3;
f = mse.value;
g = real(mse.backprop(1));
end
function [f, g] = singleWndOptSplitLFHF(weights, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
correctionWndHF = x(halfLen + halfLen + 1 : end);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S1(:, i);
    correctedTime2 = STime1(:, i) .* correctionWndHF;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse = sum(meaSqr .* normalizationVector);
f = mse.value;
g = real(mse.backprop(1));
end
function [f, g] = dualWndSplitLFMF(weights, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingMF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingMF(halfLen+1:frameSize) = wndCorrectionWeightingMF(halfLen-1:-1:2);
correctionWndTransition = x(halfLen + halfLen + 1 : end);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S1(:, i);
    correctedTime2 = STime1(:, i) .* correctionWndHF;
    correctedTime3 = STime1(:, i) .* correctionWndTransition;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    getbackCorrectedToSpectrum3 = fft(correctedTime3);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum3 .* wndCorrectionWeightingMF + getbackCorrectedToSpectrum2 .* (1 - wndCorrectionWeightingLF);
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse = sum(meaSqr .* normalizationVector);
f = mse.value;
g = real(mse.backprop(1));
end
function [f, g] = dualWndSplitLFMFHF(weights, correctionWndHF, nframes, halfLen, frameSize, hop, wnd, normalizationVector, S1, STime1, bandlimitChirp)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingMF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingMF(halfLen+1:frameSize) = wndCorrectionWeightingMF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + halfLen + 1 : halfLen + halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
correctionWndTransition = x(halfLen + halfLen + halfLen + 1 : end);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S1(:, i);
    correctedTime2 = STime1(:, i) .* correctionWndHF;
    correctedTime3 = STime1(:, i) .* correctionWndTransition;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    getbackCorrectedToSpectrum3 = fft(correctedTime3);
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum3 .* wndCorrectionWeightingMF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
    % Reflect spectrum and conjugate
    yInvShifted = ifftshift(ifft(q_fft_frameinv2)) .* wnd;
    % Overlap add
    if i > 1
        accFrame2 = accFrame2 + yInvShifted;
    else
        accFrame2 = yInvShifted;
    end
    myOut = accFrame2(1 : hop);
    accFrame2 = accFrame2(hop + 1 : frameSize);
    accFrame2(frameSize - hop + 1 : frameSize) = zeros(hop, 1);
    if i > 1
        y3(1+(i-1)*hop : i*hop) = myOut;
    else
        y3 = myOut;
    end
end
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse = sum(meaSqr .* normalizationVector);
f = mse.value;
g = real(mse.backprop(1));
end