function [coeff, f] = ltv_precomute2(fftLen, hop, fs, oct, order, HFSamplingLimit, gaussSigmaLimiting, cornerSmoothing, reqSynthesisWnd)
addpath('minFunc_2012/autoDif')
addpath('minFunc_2012/minFunc')
addpath('minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 1024;
    hop = 128;
    oct = 64;
    HFSamplingLimit = 0.7;
    order = 2;
    reqSynthesisWnd = 1;
    gaussSigmaLimiting = 1;
    cornerSmoothing = 1;
end
zp = 1;
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
ovp = fftLen / hop;
% number of points of pre and post padding used to set initial conditions
prepad = 10;
pospad = 42;
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
sigmas = (thetas1 ./ fftLen) ./ oct / pi * fftLen;
if order == 2
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
        % Eliminate oscillation around corner
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
disp('First constant Q frequnecy that undersample ' + string(f(firstUndersampling)) + ' Hz')
%% Compute reassignment frequency domain sample shifter
phaseShifter1 = exp(1i * (2*pi*(0:halfLen-1)/fftLen).');
phaseShifter2 = exp(-1i * (2*pi/fftLen).');
phaseShifter3 = exp(1i * (2*pi/fftLen).');
phaseShifter4 = exp(-1i * (2*pi*(0:halfLen-1)/fftLen).');
phaseShifter5 = exp(1i * (2*pi/fftLen).');
phaseShifter6 = exp(-1i * (2*pi/fftLen).');
%% Obtain filterbank overall response
tmp = zeros(fftLen, ovp - 1);
corrF = [];
for i = 1 : ovp - 1
    stepSize = mod(fftLen - hop * i - fftLen / 2 - 1, fftLen);
    phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
    q_fft_frameinv = ltv_1st_2ndNoSq2(phaseShifter, b, a, c1, c2, prepad, pospad, corrF, halfLen);
    %% Inverse transform
    % Virtually multiplying Hann window in time domain on frequency domain
    if reqSynthesisWnd
        q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
    end
    % Reflect spectrum and conjugate
    q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    yInvShifted = ifftshift(ifft(q_fft_frameinv)).';
    tmp(:, i) = yInvShifted;
end
systemImpulse = overlapAdd(tmp, hop);
systemImpulse = systemImpulse(1 : fftLen - hop);
if hop ~= fftLen / 2
    truncatedSystemImpulse = [systemImpulse((fftLen - hop) - (fftLen/2) : end)];
else
    truncatedSystemImpulse = [0; systemImpulse];
end
truncatedSystemImpulse = [truncatedSystemImpulse; truncatedSystemImpulse(halfLen - 1: -1 : 1)];
truncatedSystemImpulse(1) = [];
truncatedSystemImpulse = fft(truncatedSystemImpulse);
truncatedSystemImpulse = abs(truncatedSystemImpulse(1 : halfLen));
if flat
    truncatedSystemImpulse(:) = mean(truncatedSystemImpulse);
end
corrF = 1 ./ truncatedSystemImpulse;
wndCorrectionWeighting = [];
correctionWndHF = [];
theoreticalWindowShape = [];
if hop ~= 1
    %% Obtain DFT filterbank frequency response
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
    dftMtx = dftmtx(fftLen);
    dftMtx = dftMtx(1 : halfLen, :);
    cpxRes = fft( (theoreticalWindowShape .* dftMtx)' );
    cpxRes = cpxRes(1 : halfLen, :);
    dftFilterbank2 = cpxRes .* conj(cpxRes);
    overallShapeOfWeighting = mean(dftFilterbank2, 1)';
    [~, idx] = max(overallShapeOfWeighting);
    overallShapeOfWeighting(idx + 1 : end) = overallShapeOfWeighting(idx);
    overallShapeOfWeighting = overallShapeOfWeighting - min(overallShapeOfWeighting);
    if max(overallShapeOfWeighting) ~= 0
        overallShapeOfWeighting = overallShapeOfWeighting ./ max(overallShapeOfWeighting);
    end
    %% Obtain window correction
    correctionWndHF = overlapAdd(repmat(theoreticalWindowShape(end, :), ovp * 2, 1 )', hop);
    correctionWndHF = correctionWndHF(fftLen - hop + 1 : fftLen * 2 - hop);
    % The first window in constant-Q scale is widest
    if sum(theoreticalWindowShape(1, :) ./ max(theoreticalWindowShape(1, :))) < (fftLen / 5)
        error('Heavily undersampling');
    end
    correctionWndHF = 1 ./ correctionWndHF(:);
    %% Plot result after overlap-add
    if order == 2
        wndCorrectionWeighting = max(b(prepad + 1 : end - pospad + 1, :)) - b(prepad + 1 : end - pospad + 1, :);
    else
        wndCorrectionWeighting = max(a(prepad + 1 : end - pospad + 1, :)) - a(prepad + 1 : end - pospad + 1, :);
    end
    wndCorrectionWeighting = wndCorrectionWeighting - min(wndCorrectionWeighting);
    if max(wndCorrectionWeighting) ~= 0
        wndCorrectionWeighting = wndCorrectionWeighting ./ max(wndCorrectionWeighting);
    end
    fnc = @(x) fitExp(wndCorrectionWeighting, overallShapeOfWeighting, x);
    pwr = fmincon(fnc, 1, [], [], [], [], eps, Inf, [], optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'HessianFcn','objective','Display','iter'));
    wndCorrectionWeighting = 1 - (wndCorrectionWeighting .^ pwr);
    %% Iterative correction of window frequency weighting to further maximize SNR
    %% Generate chirp
    D = fftLen / fs;            % duration in seconds
    t = (0:fftLen-1) / fs;         % discrete-time axis (sec)
    f1 = 1 * fs / fftLen;
    f2 = halfLen * fs / fftLen;
    bandlimitChirp = chirp2(t, D, f1, f2)';
    %% Heuristics iterative correction of window frequency weighting to further maximize SNR
    maxIter = 100;
    bandlimitChirpPadded = [zeros(fftLen - hop, 1); bandlimitChirp; zeros(fftLen / 2, 1)];
    nframes = ovp * 2 - 1;
    % zero padding at the end to complete the last frame
    paddedLen = hop * nframes;
    bandlimitChirpPadded = [bandlimitChirpPadded; zeros(length(bandlimitChirpPadded) - paddedLen + hop * 3, 1)];
    frmIdx = 1 + (0 : nframes - 1) * hop;
    sht = fix(hop / 2);
    previousError = Inf;
    wndCorrectionWeightingPrev = wndCorrectionWeighting;
    kDeltaTFR = zeros(halfLen, ovp - 1);
    kDeltaTFRHFWnd = zeros(halfLen, ovp - 1);
    chirpTFR = zeros(halfLen, (ovp - 1) * 2);
    chirpTFRHFWnd = zeros(halfLen, (ovp - 1) * 2);
    for i = 1 : ovp * 2 - 1
        if i <= ovp - 1
            stepSize = mod((fftLen - hop * i - fftLen / 2 - 1) - sht, fftLen);
            phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
            q_fft_frameinv = ltv_1st_2ndNoSq2(phaseShifter, b, a, c1, c2, prepad, pospad, corrF, halfLen);
            if reqSynthesisWnd
                q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
            end
            kDeltaTFR(:, i) = q_fft_frameinv;
            q_fft_frameinv(halfLen+1:fftLen) = conj(q_fft_frameinv(halfLen-1:-1:2));
            correctedTime = ifft(q_fft_frameinv);
            correctedTime2 = correctedTime .* correctionWndHF;
            getbackCorrectedToSpectrum1 = fft(correctedTime2);
            kDeltaTFRHFWnd(:, i) = getbackCorrectedToSpectrum1(1 : halfLen);
        end
        frame = bandlimitChirpPadded(frmIdx(i) : frmIdx(i) + fftLen - 1);
        X = fft(fftshift(frame));
        chirpSpectrum = ltv_1st_2ndNoSq2(X(1 : halfLen), b, a, c1, c2, prepad, pospad, corrF, halfLen);
        if reqSynthesisWnd
            chirpSpectrum = 0.25 * (2 * chirpSpectrum + [conj(chirpSpectrum(2)); chirpSpectrum(1 : end - 1)] + [chirpSpectrum(2 : end); conj(chirpSpectrum(end - 1))]);
        end
        chirpTFR(:, i) = chirpSpectrum;
        chirpSpectrum(halfLen+1:fftLen) = conj(chirpSpectrum(halfLen-1:-1:2));
        correctedTime = ifft(chirpSpectrum);
        correctedTime2 = correctedTime .* correctionWndHF;
        getbackCorrectedToSpectrum2 = fft(correctedTime2);
        chirpTFRHFWnd(:, i) = getbackCorrectedToSpectrum2(1 : halfLen);
    end
    for it = 1 : maxIter
        % Overlap the approximation of system impulse
        q_fft_frameinv = kDeltaTFR .* wndCorrectionWeighting + kDeltaTFRHFWnd .* (1 - wndCorrectionWeighting);
        q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
        tmp = circshift(ifft(q_fft_frameinv), fftLen / 2);
        systemImpulse = circshift(overlapAdd(tmp, hop), sht);
        systemImpulse(1 : sht) = 0;
        systemImpulse = systemImpulse(1 : fftLen - hop);
        if hop ~= fftLen / 2
            truncatedSystemImpulse = [systemImpulse((fftLen - hop) - (fftLen/2) : end)];
        else
            truncatedSystemImpulse = [0; systemImpulse];
        end
        % Overlap the approximation of chirp signal
        chirpSpectrum = chirpTFR .* wndCorrectionWeighting + chirpTFRHFWnd .* (1 - wndCorrectionWeighting);
        % Reflect spectrum and conjugate
        chirpSpectrum(halfLen+1:fftLen,:) = conj(chirpSpectrum(halfLen-1:-1:2,:));
        tmp = circshift(ifft(chirpSpectrum), fftLen / 2);
        chirpRec = overlapAdd(tmp, hop);
        % Enforce symmetry at centre
        truncatedSystemImpulse = [truncatedSystemImpulse; truncatedSystemImpulse(halfLen - 1: -1 : 1)];
        % SAE of chirp signal
        err = chirpRec(fftLen - hop + 1 : fftLen * 2 - hop) - bandlimitChirp;
        errorCurr = sum(abs(err));
        if errorCurr >= previousError
            wndCorrectionWeighting = wndCorrectionWeightingPrev;
            break;
        else
            previousError = errorCurr;
            % disp('iteration = ' + string(it) + ', error = ' + sprintf('%1.14f', errorCurr))
        end
        truncatedSystemImpulse(1) = [];
        truncatedSystemImpulse = fft(truncatedSystemImpulse);
        truncatedSystemImpulse = abs(truncatedSystemImpulse(1 : halfLen));
        truncatedSystemImpulse = truncatedSystemImpulse - 1;
        truncatedSystemImpulse(truncatedSystemImpulse < 0) = 0;
        truncatedSystemImpulse(1) = 0;
        truncatedSystemImpulse(firstUndersampling * 2 + 1 : end) = 0;
        wndCorrectionWeightingPrev = wndCorrectionWeighting;
        wndCorrectionWeighting = wndCorrectionWeighting + truncatedSystemImpulse;
        wndCorrectionWeighting(wndCorrectionWeighting < 0) = 0;
    end
    %%
%     wndCorrectionWeightingLF = wndCorrectionWeighting;
%     wndCorrectionWeightingHF = 1 - wndCorrectionWeightingLF;
    %% Gradient based optimization
    ny = fftLen - hop + fftLen + fftLen / 2;
    halfLen = fftLen / 2 + 1;
    shiftedChirp1 = circshift(bandlimitChirp, fix(hop / 2));
    if ~isempty(firstUndersampling)
        bandlimitChirp2 = circshift(chirp2(t, D, f(2), (firstUndersampling * 2) * fs / fftLen), -fix(hop / 4))';
    else
        bandlimitChirp2 = circshift(chirp2(t, D, f(2), fs / 2 - f(2)), -fix(hop / 4))';
    end
    coeff.prepad = prepad;
    coeff.pospad = pospad;
    coeff.b = b;
    coeff.a = a;
    coeff.c1 = c1;
    coeff.c2 = c2;
    coeff.phaseShifter1 = phaseShifter1;
    coeff.phaseShifter2 = phaseShifter2;
    coeff.phaseShifter3 = phaseShifter3;
    coeff.corrF = corrF;
    coeff.fftLen = fftLen;
    coeff.halfLen = halfLen;
    coeff.hop = hop;
    coeff.fs = fs;
    [spec1, t_q] = ltv_spectrogram2(bandlimitChirp, coeff);
    [spec2, t_q] = ltv_spectrogram2(shiftedChirp1, coeff);
    [spec3, t_q] = ltv_spectrogram2(bandlimitChirp2, coeff);
    S1 = spec1(:, :);
    S2 = spec2(:, :);
    S3 = spec3(:, :);
    if reqSynthesisWnd
        S1 = 0.25 * (2 * S1 + [conj(S1(2, :)); S1(1 : end - 1, :)] + [S1(2 : end, :); conj(S1(end - 1, :))]);
        S2 = 0.25 * (2 * S2 + [conj(S2(2, :)); S2(1 : end - 1, :)] + [S2(2 : end, :); conj(S2(end - 1, :))]);
        S3 = 0.25 * (2 * S3 + [conj(S3(2, :)); S3(1 : end - 1, :)] + [S3(2 : end, :); conj(S3(end - 1, :))]);
    end
    S1(halfLen+1:fftLen,:) = conj(S1(halfLen-1:-1:2,:));
    S2(halfLen+1:fftLen,:) = conj(S2(halfLen-1:-1:2,:));
    S3(halfLen+1:fftLen,:) = conj(S3(halfLen-1:-1:2,:));
    STime1 = ifft(S1);
    STime2 = ifft(S2);
    STime3 = ifft(S3);
    normalizationVector = ones(fftLen, 1) / fftLen;
    wndCorrectionWeightingLF = wndCorrectionWeighting;
    wndCorrectionWeightingHF = 1 - wndCorrectionWeightingLF;
    weights1 = [wndCorrectionWeightingLF; wndCorrectionWeightingHF];
    getbackCorrectedToSpectrum1 = fft(STime1 .* correctionWndHF);
    getbackCorrectedToSpectrum2 = fft(STime2 .* correctionWndHF);
    getbackCorrectedToSpectrum3 = fft(STime3 .* correctionWndHF);
    func1 = @(x) singleWndSplitLFHF(x, halfLen, fftLen, hop, normalizationVector, S1, getbackCorrectedToSpectrum1, bandlimitChirp, S2, getbackCorrectedToSpectrum2, shiftedChirp1, S3, getbackCorrectedToSpectrum3, bandlimitChirp2);
    opt.optTol = 1e-14;
    opt.progTol = 1e-14;
    opt.MaxFunEvals = 1000;
    opt.MaxIter = 1000;
    opt.Method = 'pnewton0';
    weightsOpt = minFunc(func1, weights1, opt);
    coeff.wndCorrectionWeightingLF = movmedian(weightsOpt(1 : halfLen), 10);
    coeff.wndCorrectionWeightingHF = movmedian(weightsOpt(halfLen + 1 : halfLen + halfLen), 10);

%     wndCorrectionWeightingLF = coeff.wndCorrectionWeightingLF;
%     wndCorrectionWeightingHF = coeff.wndCorrectionWeightingHF;
%     testedWndowShape = zeros(size(theoreticalWindowShape));
%     correctionMatrix = zeros(halfLen, halfLen);
%     X = zeros(halfLen, 1);
%     for idx = 1 : halfLen
%         X(:) = 0;
%         X(idx) = 1;
%         q_fft_frameinv = ltv_1st_2ndNoSq2(X, b, a, c1, c2, prepad, pospad, corrF, halfLen);
%         q_fft_frameinv2 = q_fft_frameinv;
%         q_fft_frameinv2(halfLen+1:fftLen,:) = conj(q_fft_frameinv2(halfLen-1:-1:2,:));
%         correctedTime = ifft(q_fft_frameinv2);
%         correctedTime2 = correctedTime .* correctionWndHF;
%         getbackCorrectedToSpectrum2 = fft(correctedTime2);
%         % Frequency weighted interpolation between original and windowed spectra
%         q_fft_frameinv = q_fft_frameinv .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2(1 : halfLen) .* wndCorrectionWeightingHF;
%         if reqSynthesisWnd
%             % Virtually multiplying Hann window in time domain on frequency domain
%             q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
%         end
%         wndCentred = circshift(q_fft_frameinv, -idx + fftLen / 4, 1);
%         gTime = ifftshift(abs(fft(wndCentred, fftLen)))';
%         correctionWnd = overlapAdd(repmat(gTime, ovp * 2, 1 )', hop);
%         correctionWnd = correctionWnd(fftLen - hop + 1 : fftLen * 2 - hop);
%         halfCorrWnd = 1 ./ correctionWnd(1 : halfLen);
%         halfCorrWnd(halfLen+1:fftLen) = conj(halfCorrWnd(halfLen-1:-1:2));
%         %     halfCorrWnd = correctionWndHF;
%         invHalfCorr = ifft(halfCorrWnd);
%         invHalfCorr = ifftshift(invHalfCorr);
%         g2 = circshift(invHalfCorr, idx - fftLen / 2 - 1);
%         g2 = g2(1 : halfLen);
%         testedWndowShape(idx, :) = gTime;
%         correctionMatrix(idx, :) = g2;
%     end
    correctionMatrix = [];
else
    correctionMatrix = [];
end
%% Compute correction for synchrosqueezing transform
corrS = ones(halfLen, 1);
max1 = -Inf;
max2 = -Inf;
cnt = 0;
firstBestSol = zeros(halfLen, 1);
it = 0;
if hop == 1
    target = zeros(fftLen, 1);
    target(fftLen / 2) = 1;
    kDeltaSpectrum = zeros(halfLen, fftLen/hop);
    for i=1:fftLen/hop
        stepSize = mod(fftLen - hop * i - fftLen / 2, fftLen);
        phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
        S2 = ltv_1st_2nd(phaseShifter, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, corrS, 1);
        kDeltaSpectrum(:, i) = S2(:, 2);
    end
    while (1)
        S2 = kDeltaSpectrum .* corrS;
        synsq = sum(S2(:, :), 1) + conj( sum(S2(2 : fftLen / 2, :), 1) );
        synsq = real(synsq(:));
        synsqSpec = fft(synsq);
        synsqSpec = abs(synsqSpec(1 : halfLen));
        SNR = 10*log10(sum(abs(target).^2)/sum(abs(target-synsq).^2)); % Calculate the SNR
        it = it + 1;
        if it > 1000
            break;
        end
        if SNR > max1
            max2 = max1;
            max1 = SNR;
            firstBestSol(:) = corrS;
            cnt = 0;
        elseif SNR > max2 && SNR < max1
            max2 = SNR;
            cnt = 0;
        else
            % If current solution worse than before, we still run for
            % hundred rounds to capture
            cnt = cnt + 1;
            if cnt > 100 % If solution quality keep decreasing for more than hundred round, then terminate
                corrS(:) = firstBestSol;
                break;
            end
        end
        corrS = corrS .* 1 ./ synsqSpec;
    end
    %      fvtool(synsq);
end
%% Save coefficient
coeff.correctionMatrix = correctionMatrix;
coeff.prepad = prepad;
coeff.pospad = pospad;
coeff.b = b;
coeff.a = a;
coeff.c1 = c1;
coeff.c2 = c2;
coeff.phaseShifter1 = phaseShifter1;
coeff.phaseShifter2 = phaseShifter2;
coeff.phaseShifter3 = phaseShifter3;
coeff.corrF = corrF;
coeff.corrS = corrS;
coeff.correctionWnd = correctionWndHF;
coeff.wndCorrectionWeighting = wndCorrectionWeighting;
coeff.fftLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = hop;
coeff.fs = fs;
coeff.theoreticalWindowShape = theoreticalWindowShape;
end
function [f, g, h] = fitExp(corrF, cdf, raise)
f = mean( ( (corrF .^ raise) - cdf) .^ 2, 1); % Cost
if nargout > 1
    J1 = -(2*corrF.^raise.*log(corrF) .* (cdf - corrF.^raise)); % Gradient
    J1(isnan(J1)) = 0;
    g = mean(J1, 1);
end
if nargout > 2
    H1 = 4 * [corrF.^(2*raise).*log(corrF).^2; -corrF.^raise.*log(corrF).^2.*(cdf - corrF.^raise)]; % Hessian
    H1(isnan(H1)) = 0;
    h = mean(H1, 1);
end
end
function x = chirp2(t,t1,f0,f1)
beta = (f1-f0)./t1;
x = cos(2*pi * ( 0.5* beta .* (t .* t) + f0 * t));
end
function [f, g] = singleWndSplitLFHF(weights, halfLen, frameSize, hop, normalizationVector, S1, getbackCorrectedToSpectrum1, bandlimitChirp, S2, getbackCorrectedToSpectrum2, shiftedChirp1, S3, getbackCorrectedToSpectrum3, shiftedChirp2)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
%% Inverse transform
% Frequency weighted interpolation between original and windowed spectra
q_fft_frameinv2 = S1 .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum1 .* wndCorrectionWeightingHF;
mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
y3 = fold(mtx, frameSize, hop);
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse1 = sum(meaSqr .* normalizationVector);
%% Inverse transform
% Frequency weighted interpolation between original and windowed spectra
q_fft_frameinv2 = S2 .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
y3 = fold(mtx, frameSize, hop);
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - shiftedChirp1) .^ 2;
mse2 = sum(meaSqr .* normalizationVector);
%% Inverse transform
% Frequency weighted interpolation between original and windowed spectra
q_fft_frameinv2 = S3 .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum3 .* wndCorrectionWeightingHF;
mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
y3 = fold(mtx, frameSize, hop);
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - shiftedChirp2) .^ 2;
mse3 = sum(meaSqr .* normalizationVector);
mse = (mse1 + mse2 + mse3) / 3;
f = mse.value;
g = real(mse.backprop(1));
end
function [f, g] = dualWndSplitLFMFHF(weights, correctionWndHF, correctionWndTransition, halfLen, frameSize, hop, normalizationVector, S1, STime1, bandlimitChirp, S2, STime2, shiftedChirp1, S3, STime3, shiftedChirp2, optimizableMFWnd)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingMF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingMF(halfLen+1:frameSize) = wndCorrectionWeightingMF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + halfLen + 1 : halfLen + halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
if optimizableMFWnd
    correctionWndTransition = x(halfLen + halfLen + halfLen + 1 : end);
end
%% Inverse transform
correctedTime2 = STime1 .* correctionWndHF;
correctedTime3 = STime1 .* correctionWndTransition;
getbackCorrectedToSpectrum2 = fft(correctedTime2);
getbackCorrectedToSpectrum3 = fft(correctedTime3);
% Frequency weighted interpolation between original and windowed spectra
q_fft_frameinv2 = S1 .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum3 .* wndCorrectionWeightingMF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
y3 = fold(mtx, frameSize, hop);
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - bandlimitChirp) .^ 2;
mse1 = sum(meaSqr .* normalizationVector);
%% Inverse transform
correctedTime2 = STime2 .* correctionWndHF;
correctedTime3 = STime2 .* correctionWndTransition;
getbackCorrectedToSpectrum2 = fft(correctedTime2);
getbackCorrectedToSpectrum3 = fft(correctedTime3);
% Frequency weighted interpolation between original and windowed spectra
q_fft_frameinv2 = S2 .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum3 .* wndCorrectionWeightingMF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
y3 = fold(mtx, frameSize, hop);
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - shiftedChirp1) .^ 2;
mse2 = sum(meaSqr .* normalizationVector);
%% Inverse transform
correctedTime2 = STime3 .* correctionWndHF;
correctedTime3 = STime3 .* correctionWndTransition;
getbackCorrectedToSpectrum2 = fft(correctedTime2);
getbackCorrectedToSpectrum3 = fft(correctedTime3);
% Frequency weighted interpolation between original and windowed spectra
q_fft_frameinv2 = S3 .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum3 .* wndCorrectionWeightingMF + getbackCorrectedToSpectrum2 .* wndCorrectionWeightingHF;
mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
y3 = fold(mtx, frameSize, hop);
y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
meaSqr = (y - shiftedChirp2) .^ 2;
mse3 = sum(meaSqr .* normalizationVector);
mse = (mse1 + mse2 + mse3) / 3;
f = mse.value;
g = real(mse.backprop(1));
end
function [S, t] = ltv_spectrogram2(x, coeff)
x = x(:);
x = [zeros(coeff.fftLen - coeff.hop, 1); x; zeros(coeff.fftLen / 2, 1)];
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.fftLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
x = [x; zeros(length(x) - paddedLen + coeff.hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
S = zeros(nframes, coeff.fftLen/2+1, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.fftLen - 1);
    X = fft(fftshift(frame));
    S(i, :) = ltv_1st_2ndNoSq2(X(1 : coeff.halfLen), coeff.b, coeff.a, coeff.c1, coeff.c2, coeff.prepad, coeff.pospad, coeff.corrF, coeff.halfLen);
end
t = (frmIdx - 1 + coeff.fftLen / 2) / coeff.fs;
S = permute(S, [2, 1, 3]);
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2, prepad, pospad, corrF, halfLen)
%% Hann in frequency domain with input being shifted by fftshift
specHann = 2 * dftSpec + [conj(dftSpec(2)); dftSpec(1 : end - 1)] + [dftSpec(2 : end); conj(dftSpec(end - 1))];
%% Hann in frequency domain with 1 sample delayed input being shifted by fftshift
% Hann windowed
x_fft1 = [conj(specHann(prepad + 1 : -1 : 2, :)); specHann; conj(specHann(halfLen - 1 : -1 : (halfLen - pospad + 1), :))] / 4;
% Rectangular windowed
% x_fft1 = [conj(dftSpec(prepad + 1 : -1 : 2, :)); dftSpec; conj(dftSpec(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
%% Gaussian windowing
tmp = zeros(size(x_fft1, 1), 1, 'like', x_fft1);
q_fft_frame = ltv1Slice(x_fft1, tmp, b, a, c1, c2);
% Remove periodic padding
spec = q_fft_frame(prepad + 1 : end - pospad + 1);
if ~isempty(corrF)
    spec = spec .* corrF;
end
spec(1, :) = real(spec(1, :));
spec(end, :) = real(spec(end, :));
end