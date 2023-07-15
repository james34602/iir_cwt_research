function [coeff, f] = ltv_precomute2(fftLen, hop, fs, oct, order, HFSamplingLimit, gaussSigmaLimiting, cornerSmoothing, reqSynthesisWnd, sparCutOff)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
addpath('gradients/')
if nargin < 1
    fs = 48000;
    fftLen = 4096;
    hop = 512;
    oct = 64;
    HFSamplingLimit = 1.5;
    order = 2;
    reqSynthesisWnd = 1;
    gaussSigmaLimiting = 1;
    cornerSmoothing = 1;
    sparCutOff = 1e-7;
end
zp = 1;
halfLen = getFFTHalfLen(fftLen);
ovp = fftLen / hop;
% number of points of pre and post padding used to set initial conditions
prepad = fix(fftLen / 102.4);
pospad = fix(fftLen / 102.4) + 2;
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
analysisWnd = tukeywin_periodic(fftLen, 0.5);
analysisWnd = hann(fftLen, 'periodic') .^ (1 / reqSynthesisWnd);
synWnd = hann(fftLen, 'periodic');
chopedWnd1 = analysisWnd(halfLen : end);
chopedWnd2 = synWnd(halfLen : end);
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
hHopPt = hHopPtCplx .* conj(hHopPtCplx) * chopedWnd1(div);
hHopPt = hHopPt * (chopedWnd2(div) .^ reqSynthesisWnd);
top = max(hHopPt); down = min(hHopPt);
firstUndersampling = find(hHopPt(prepad + 1 : end - pospad + 1) <= ((top + down) / 2), 1, 'first');
% firstUndersampling = find((hHopPt(prepad + 1 : end - pospad + 1) * 2) <= 1, 1, 'first');
firstUndersampling = min(firstUndersampling, fix(fftLen / 4 - 1) - prepad);
if ~isempty(prepad + firstUndersampling - 1)
    thetas1(prepad + firstUndersampling - 1 : halfLen) = thetas1(prepad + firstUndersampling - 1);
    thetas1(1 : prepad) = thetas1(prepad * 2 + 1 : -1 : prepad + 2);
    thetas1(halfLen + 1 : end) = thetas1(halfLen-1:-1:halfLen - prepad - pospad + 1);
    if cornerSmoothing
        % Eliminate oscillation around corner
        time = 0.026 * fftLen / zp; % More elements in array less smoothing is needed
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
%% Obtain filterbank overall response
tmp = zeros(fftLen, ovp - 1);
corrF = [];
imp = zeros(fftLen, 1);
for i = 1 : ovp - 1
    stepSize = mod(fftLen - hop * i, fftLen);
    imp(:) = 0;
    imp(stepSize) = 1;
    phaseShifter = fft(fftshift(imp .* analysisWnd));
    phaseShifter = phaseShifter(1 : halfLen);
    q_fft_frameinv = ltv_1st_2ndNoSq2(phaseShifter, b, a, c1, c2, prepad, pospad, corrF, halfLen);
    %% Inverse transform
    % Virtually multiplying Hann window in time domain on frequency domain
    for nWnd = 1 : reqSynthesisWnd
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
%% Obtain DFT filterbank frequency response
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd1.';
h2 = h2 .* (chopedWnd2 .^ reqSynthesisWnd).';
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
    warning('Heavily undersampling');
end
correctionWndHF = 1 ./ correctionWndHF(:);
% [h,w] = freqz(correctionWndHF(1 : halfLen),1,halfLen);
% figure(8)
% plot(real(fft(correctionWndHF)));
% hold on
% plot(imag(fft(correctionWndHF)))
% hold off
% axis tight
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
        stepSize = mod(fftLen - hop * i - sht, fftLen);
        imp(:) = 0;
        imp(stepSize) = 1;
        phaseShifter = fft(fftshift(imp .* analysisWnd));
        phaseShifter = phaseShifter(1 : halfLen);
        q_fft_frameinv = ltv_1st_2ndNoSq2(phaseShifter, b, a, c1, c2, prepad, pospad, corrF, halfLen);
        for nWnd = 1 : reqSynthesisWnd
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
    X = fft(fftshift(frame .* analysisWnd));
    chirpSpectrum = ltv_1st_2ndNoSq2(X(1 : halfLen), b, a, c1, c2, prepad, pospad, corrF, halfLen);
    for nWnd = 1 : reqSynthesisWnd
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
coeff.analysisWnd = analysisWnd;
coeff.prepad = prepad;
coeff.pospad = pospad;
coeff.b = b;
coeff.a = a;
coeff.c1 = c1;
coeff.c2 = c2;
coeff.phaseShifter1 = phaseShifter1;
coeff.corrF = corrF;
coeff.fftLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = hop;
coeff.fs = fs;
numDemoSignal = 3;
sigs = cell(numDemoSignal, 1);
wndCorrectionWeightingLF = wndCorrectionWeighting;
wndCorrectionWeightingHF = 1 - wndCorrectionWeightingLF;
weights1 = [wndCorrectionWeightingLF; wndCorrectionWeightingHF];
for idx = 1 : numDemoSignal
    if idx == 1
        target = bandlimitChirp;
    elseif idx == 2
        target = circshift(bandlimitChirp, fix(hop / 2));
    else
        if ~isempty(firstUndersampling)
            target = circshift(chirp2(t, D, f(2), (firstUndersampling * 2) * fs / fftLen), -fix(hop / 4))';
        else
            target = circshift(chirp2(t, D, f(2), fs / 2 - f(2)), -fix(hop / 4))';
        end
    end
    [spec1, t_q] = ltv_spectrogram2(target, coeff);
    S1 = spec1(:, :);
    for nWnd = 1 : reqSynthesisWnd
        S1 = 0.25 * (2 * S1 + [conj(S1(2, :)); S1(1 : end - 1, :)] + [S1(2 : end, :); conj(S1(end - 1, :))]);
    end
    S1(halfLen+1:fftLen,:) = conj(S1(halfLen-1:-1:2,:));
    STime1 = ifft(S1);
    getbackCorrectedToSpectrum1 = fft(STime1 .* correctionWndHF);
    pack.S = S1;
    pack.SRe = real(S1(1:halfLen, :));
    pack.SIm = imag(S1(1:halfLen, :));
    pack.getbackCorrectedToSpectrum = getbackCorrectedToSpectrum1;
    pack.getbackCorrectedToSpectrumRe = real(getbackCorrectedToSpectrum1(1:halfLen, :));
    pack.getbackCorrectedToSpectrumIm = imag(getbackCorrectedToSpectrum1(1:halfLen, :));
    pack.target = target;
    if idx == 1
        sg.SRe = zeros([size(pack.SRe), numDemoSignal]);
        sg.SIm = zeros([size(pack.SIm), numDemoSignal]);
        sg.getbackCorrectedToSpectrumRe = zeros([size(pack.getbackCorrectedToSpectrumRe), numDemoSignal]);
        sg.getbackCorrectedToSpectrumIm = zeros([size(pack.getbackCorrectedToSpectrumIm), numDemoSignal]);
        sg.target = zeros([size(target), numDemoSignal]);
        sg.target = squeeze(sg.target);
    end
    sg.SRe(:, :, idx) = pack.SRe;
    sg.SIm(:, :, idx) = pack.SIm;
    sg.getbackCorrectedToSpectrumRe(:, :, idx) = pack.getbackCorrectedToSpectrumRe;
    sg.getbackCorrectedToSpectrumIm(:, :, idx) = pack.getbackCorrectedToSpectrumIm;
    sg.target(:, idx) = target;
    sigs{idx, 1} = pack;
end
func1 = @(x) singleWndSplitLFHF(x, halfLen, fftLen, hop, sigs);
func2 = @(x) singleWndSplitLFHFPaged(x, halfLen, fftLen, hop, sg);
func3 = @(x) singleWndSplitLFHFEnergyConservation(x, halfLen, fftLen, hop, sigs);
opt.optTol = 1e-14;
opt.progTol = 1e-14;
opt.MaxIter = 1000;
opt.MaxFunEvals = opt.MaxIter * 10;
opt.Method = 'pnewton0';
dbg = 1;
if dbg
    weightsOpt = minFunc(func2, weights1, opt);
    coeff.wndCorrectionWeightingLF = movmedian(weightsOpt(1 : halfLen), 10);
    coeff.wndCorrectionWeightingHF = movmedian(weightsOpt(halfLen + 1 : halfLen + halfLen), 10);
else
    load('matlab.mat');
end
%%
wndCorrectionWeightingLF = coeff.wndCorrectionWeightingLF;
wndCorrectionWeightingHF = coeff.wndCorrectionWeightingHF;
%% DFT and IDFT matrix
% Wdft = fft(eye(fftLen));
% WdftRe = real(Wdft(1 : halfLen, :));
% WdftIm = imag(Wdft(1 : halfLen, :));
% Winvdft = conj(Wdft) / fftLen;
% WinvdftRe = real(Winvdft);
% WinvdftIm = imag(Winvdft); WinvdftIm(:, halfLen + 1 : end) = -WinvdftIm(:, halfLen + 1: end);
%% Optimization parameters
opt.optTol = eps;
opt.progTol = eps;
opt.MaxIter = 1000;
opt.MaxFunEvals = opt.MaxIter * 2;
opt.Method = 'pnewton0';
%% Prepare signals
desiredFrameCount = fftLen / 16;
cpLen = ceil((desiredFrameCount + (fftLen / hop / 2)) * hop - fftLen - hop - fftLen / 2);
D = cpLen / fs;t = (0:cpLen-1) / fs;f1 = 1 * fs / fftLen;f2 = halfLen * fs / fftLen;longChirp = chirp2(t, D, f1, f2)';
numDemoSignal = 16;
sigs = cell(numDemoSignal, 1);
pack = [];
for idx = 1 : numDemoSignal
    target = circshift(longChirp, fix((idx - 1) * (hop / numDemoSignal)));
    [S1, t_q] = ltv_spectrogram2(target, coeff);
    for nWnd = 1 : reqSynthesisWnd
        S1 = 0.25 * (2 * S1 + [conj(S1(2, :)); S1(1 : end - 1, :)] + [S1(2 : end, :); conj(S1(end - 1, :))]);
    end
    S1(halfLen+1:fftLen,:) = conj(S1(halfLen-1:-1:2,:));
    STime1 = ifft(S1);
    getbackCorrectedToSpectrum1 = fft(STime1 .* correctionWndHF);
    S1 = S1(1 : halfLen, :) .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum1(1 : halfLen, :) .* wndCorrectionWeightingHF;
    % Optimize real matrix with extended spectral(Extended spectral simulate circular convolution)
    ext = [conj(S1(prepad + 1 : -1 : 2, :)); S1; conj(S1(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
    if idx == 1
        sg.SRe = zeros([size(ext), numDemoSignal]);
        sg.SIm = zeros([size(ext), numDemoSignal]);
        sg.target = zeros([size(target), numDemoSignal]);
        sg.target = squeeze(sg.target);
    end
    sg.SRe(:, :, idx) = real(ext);
    sg.SIm(:, :, idx) = imag(ext);
    sg.target(:, idx) = target;
%     pack.SRe = sg.SRe(:, :, idx);
%     pack.SIm = sg.SIm(:, :, idx);
%     pack.target = target;
%     sigs{idx} = pack;
end
%%
identity = eye(prepad + halfLen + pospad - 1);identity = identity(prepad + 1 : end - pospad + 1, :);
identity = identity(:);
func2_ = @(x) corrRealMtxFastExtendedPagedHC(x, halfLen, fftLen, prepad, pospad, hop, sg, []);
% func2_ = @(x) corrRealMtxFastExtended(x, halfLen, fftLen, prepad, pospad, hop, sigs, []);
% [fval1_, grad1_, preview1_] = func2_(identity);
weightsOpt3 = minFunc(func2_, identity, opt);
enableSparsification = true;
enablePostOptimization = true;
%% Sparsification
gWeights = abs(weightsOpt3) > sparCutOff;
sparsity = 1 - sum(gWeights(:)) / numel(weightsOpt3)
if enableSparsification
    weightsOpt3 = weightsOpt3 .* gWeights;
    func2_ = @(x) corrRealMtxFastExtendedPagedHC(x, halfLen, fftLen, prepad, pospad, hop, sg, gWeights);
    opt.MaxIter = 100;
    opt.MaxFunEvals = opt.MaxIter * 2;
    weightsOpt3 = minFunc(func2_, weightsOpt3, opt);
end
correctionMatrix = reshape(weightsOpt3, halfLen, halfLen + prepad + pospad - 1);
%% Squeeze more juice
if enablePostOptimization
    sigs = cell(numDemoSignal, 1);
    for idx = 1 : numDemoSignal
        target = circshift(longChirp, fix((idx - 1) * (hop / numDemoSignal)));
        [spec1, t_q] = ltv_spectrogram2(target, coeff);
        S1 = spec1(:, :);
        for nWnd = 1 : reqSynthesisWnd
            S1 = 0.25 * (2 * S1 + [conj(S1(2, :)); S1(1 : end - 1, :)] + [S1(2 : end, :); conj(S1(end - 1, :))]);
        end
        S1(halfLen+1:fftLen,:) = conj(S1(halfLen-1:-1:2,:));
        STime1 = ifft(S1);
        getbackCorrectedToSpectrum1 = fft(STime1 .* correctionWndHF);
        pack.SRe = real(S1(1:halfLen, :));
        pack.SIm = imag(S1(1:halfLen, :));
        pack.SRe = [pack.SRe(prepad + 1 : -1 : 2, :); pack.SRe; pack.SRe(halfLen - 1 : -1 : (halfLen - pospad + 1), :)];
        pack.SIm = [-pack.SIm(prepad + 1 : -1 : 2, :); pack.SIm; -pack.SIm(halfLen - 1 : -1 : (halfLen - pospad + 1), :)];
        pack.getbackCorrectedToSpectrumRe = real(getbackCorrectedToSpectrum1(1:halfLen, :));
        pack.getbackCorrectedToSpectrumIm = imag(getbackCorrectedToSpectrum1(1:halfLen, :));
        pack.getbackCorrectedToSpectrumRe = [pack.getbackCorrectedToSpectrumRe(prepad + 1 : -1 : 2, :); pack.getbackCorrectedToSpectrumRe; pack.getbackCorrectedToSpectrumRe(halfLen - 1 : -1 : (halfLen - pospad + 1), :)];
        pack.getbackCorrectedToSpectrumIm = [-pack.getbackCorrectedToSpectrumIm(prepad + 1 : -1 : 2, :); pack.getbackCorrectedToSpectrumIm; -pack.getbackCorrectedToSpectrumIm(halfLen - 1 : -1 : (halfLen - pospad + 1), :)];
        pack.target = target;
        sigs{idx, 1} = pack;
    end
    func1 = @(x) singleWndSplitLFHFMtx(x, fftLen, prepad, pospad, halfLen, hop, sigs, correctionMatrix);
    wt = [coeff.wndCorrectionWeightingLF; coeff.wndCorrectionWeightingHF];
    opt.optTol = eps * 0.01;
    opt.progTol = eps * 0.01;
    opt.MaxIter = 100;
    opt.MaxFunEvals = opt.MaxIter * 2;
    opt.Method = 'pnewton0';
    weightsOpt = minFunc(func1, wt, opt);
    coeff.wndCorrectionWeightingLF = weightsOpt(1 : halfLen);
    coeff.wndCorrectionWeightingHF = weightsOpt(halfLen + 1 : halfLen + halfLen);
end
%% Optimize SVF filter with extended spectral(Extended spectral simulate circular convolution)
% boundedZero = 1;
% extendedLen = halfLen + prepad + pospad - 1;
% if boundedZero
%     b0 = zeros(extendedLen, 1);b1 = zeros(extendedLen, 1);b2 = zeros(extendedLen, 1);
%     for idx = 1 : extendedLen
%         while (1)
%             b = [1, randn(1, 2)];
%             rt1 = roots(b);
%             if any(abs(rt1) > 1)
%                 continue;
%             else
%                 break;
%             end
%         end
%         b0(idx) = b(1);b1(idx) = b(2);b2(idx) = b(3);
%     end
% else
%     b0 = randn(extendedLen, 1);b1 = randn(extendedLen, 1);b2 = randn(extendedLen, 1);
% end
% a1 = zeros(extendedLen, 1);a2 = zeros(extendedLen, 1);
% a1(:) = b1(:);a2(:) = b2(:);
% c1 = a1 + 2;
% c2 = (1 + a1 + a2) ./ c1;
% d1 = (2 * b0 + b1) ./ c1;
% d2 = (b0 + b1 + b2) ./ (c1 .* c2);
% poles = [b0(:); d1(:); d2(:); c1(:); c2(:)];
% func2_ = @(x) corrRealSVFExtended(x, halfLen, fftLen, prepad, pospad, hop, sigs, []);
% [fval2_, grad2_, preview2_] = func2_(poles);
% opt.Method = 'lbfgs'; % pcg, lbfgs
% weightsOpt2_ = minFunc(func2_, poles, opt);
% disp('')
%% Optimize complex matrix with extended spectral
% identity = [eye(prepad + halfLen + pospad - 1), zeros(prepad + halfLen + pospad - 1)]; identity = identity(:);
% func5 = @(x) corrCplxMtxFastExtended(x, halfLen, fftLen, prepad, pospad, hop, sigs, []);
% [fval5, grad5, preview5] = func5(identity);
% weightsOpt5 = minFunc(func5, identity, opt);
%% Save coefficient
coeff.correctionMatrix = sparse(correctionMatrix);
coeff.prepad = prepad;
coeff.pospad = pospad;
coeff.b = b;
coeff.a = a;
coeff.c1 = c1;
coeff.c2 = c2;
coeff.phaseShifter1 = phaseShifter1;
coeff.corrF = corrF;
coeff.correctionWnd = correctionWndHF;
coeff.fftLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = hop;
coeff.fs = fs;
coeff.reqSynthesisWnd = reqSynthesisWnd;
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
    X = fft(fftshift(frame .* coeff.analysisWnd));
    S(i, :) = ltv_1st_2ndNoSq2(X(1 : coeff.halfLen), coeff.b, coeff.a, coeff.c1, coeff.c2, coeff.prepad, coeff.pospad, coeff.corrF, coeff.halfLen);
end
t = (frmIdx - 1 + coeff.fftLen / 2) / coeff.fs;
S = permute(S, [2, 1, 3]);
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2, prepad, pospad, corrF, halfLen)
% Rectangular windowed
x_fft1 = [conj(dftSpec(prepad + 1 : -1 : 2, :)); dftSpec; conj(dftSpec(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
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
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function w = tukeywin_periodic(N, r)
if r <= 0
    w = ones(N,1);
elseif r >= 1
    w = hann(N, 'periodic');
else
    L = N + 1;
    t = linspace(0,1,L)';
    % Defines period of the taper as 1/2 period of a sine wave.
    per = r/2; 
    tl = floor(per*(L-1))+1;
    th = L-tl+1;
    % Window is defined in three sections: taper, constant, taper
    w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th : end - 1) - 1 + per)))/2)];
end
end