function [coeff, f] = ltv_precomute(fftLen, hop, fs, oct, order, HFSamplingLimit, gaussSigmaLimiting, cornerSmoothing, reqSynthesisWnd)
if nargin < 1
    fs = 48000;
    fftLen = 256;
    hop = 16;
    oct = 33;
    HFSamplingLimit = 0.7;
    order = 2;
    reqSynthesisWnd = 1;
    gaussSigmaLimiting = 1;
    cornerSmoothing = 1;
end
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
%% Obtain filterbank overall response
systemImpulseAcc = zeros(fftLen, 1);
systemImpulse = zeros(fftLen - hop, 1);
corrF = [];
corrS = [];
cnt = 1;
for i=1:fftLen/hop - 1
    stepSize = mod(fftLen - hop * i - fftLen / 2 - 1, fftLen);
    phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
    S2 = ltv_1st_2ndNoSq(phaseShifter, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, [], []);
    %% Inverse transform
    q_fft_frameinv = S2(:, 1);
    % Virtually multiplying Hann window in time domain on frequency domain
    if reqSynthesisWnd
        q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
    end
    % Reflect spectrum and conjugate
    q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    yInvShifted = ifftshift(ifft(q_fft_frameinv));
    % Overlap add
    systemImpulseAcc = systemImpulseAcc + yInvShifted;
    myOut = systemImpulseAcc(1 : hop);
    systemImpulseAcc = [systemImpulseAcc(hop + 1 : end); zeros(hop, 1)];
    systemImpulse(cnt : cnt + hop - 1) = myOut;
    cnt = cnt + hop;
end
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
wndCorrectionWeighting = [];
correctionWnd = [];
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
    % plot(theoreticalWindowShape')
    % axis tight
     dftMtx = dftmtx(fftLen);
     dftMtx = dftMtx(1 : halfLen, :);
     cpxRes = fft( (theoreticalWindowShape .* dftMtx)' );
     cpxRes = cpxRes(1 : halfLen, :);
     dftFilterbank2 = cpxRes .* conj(cpxRes);
    overallShapeOfWeighting = mean(dftFilterbank2, 1)';
%     overallShapeOfWeighting = 1 ./ mean(theoreticalWindowShape, 2);
%     if overallShapeOfWeighting(end) < overallShapeOfWeighting(1)
%         overallShapeOfWeighting = -overallShapeOfWeighting + mean(overallShapeOfWeighting);
%     end
    [~, idx] = max(overallShapeOfWeighting);
    overallShapeOfWeighting(idx + 1 : end) = overallShapeOfWeighting(idx);
    overallShapeOfWeighting = overallShapeOfWeighting - min(overallShapeOfWeighting);
    if max(overallShapeOfWeighting) ~= 0
        overallShapeOfWeighting = overallShapeOfWeighting ./ max(overallShapeOfWeighting);
    end
    %% Obtain window correction
    correctionWnd = zeros(1, size(theoreticalWindowShape, 2));
    % figure(1)
    % clf
    % hold on
    for idx = 1 : fftLen / hop
        shifted2 = circshift(theoreticalWindowShape(end, :), (idx - 1) * hop, 2);
        %     nd = shifted(firstUndersampling, :) ./ max(shifted(firstUndersampling, :));
        %     plot(nd);
        correctionWnd = correctionWnd + shifted2;
    end
    % The first window in constant-Q scale is widest
    if sum(theoreticalWindowShape(1, :) ./ max(theoreticalWindowShape(1, :))) < (fftLen / 5)
        error('Heavily undersampling');
    end
    correctionWnd = 1 ./ correctionWnd(:);
    % hold off
    %% Plot result after overlap-add
    % summ = zeros(size(theoreticalWindowShape, 1), size(theoreticalWindowShape, 2));
    % figure(1)
    % clf
    % hold on
    % for idx = 1 : fftLen / hop
    %     shifted = circshift(theoreticalWindowShape, (idx - 1) * hop, 2);
    %     nd = shifted(firstUndersampling, :) ./ max(shifted(firstUndersampling, :));
    %     plot(nd);
    %     summ = summ + shifted;
    % end
    % hold off
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
    pwr = 1;
    % [pwr, ~] = cgtrust(pwr, fnc, [eps, 0.1, 1000, 1000]);
    pwr = fminunc(fnc, pwr, optimoptions('fminunc','SpecifyObjectiveGradient',true, 'Display','off'));
    wndCorrectionWeighting = 1 - (wndCorrectionWeighting .^ pwr);
    %% Iterative correction of window frequency weighting to further maximize SNR
    maxIter = 100;
    %% Generate chirp
    D = fftLen / fs;            % duration in seconds
    t = (0:fftLen-1) / fs;         % discrete-time axis (sec)
    if ~isempty(firstUndersampling)
        f1 = fix(firstUndersampling / 3) * fs / fftLen;
        f2 = fix(firstUndersampling * 2) * fs / fftLen;
    else
        f1 = 1 * fs / fftLen;
        f2 = halfLen * fs / fftLen;
    end
    bandlimitChirp = chirp2(t, D, f1, f2)';
    bandlimitChirpPadded = [zeros(fftLen, 1); bandlimitChirp; zeros(fftLen / 2, 1)];
    ny = length(bandlimitChirpPadded);
    nframes = fix(ceil(ny/hop) - fftLen / hop / 2);
    % zero padding at the end to complete the last frame
    paddedLen = hop * nframes;
    bandlimitChirpPadded = [bandlimitChirpPadded; zeros(length(bandlimitChirpPadded) - paddedLen + hop * 2, 1)];
    frmIdx = 1 + (0 : nframes - 1) * hop;
    sht = fix(hop / 2);
    systemImpulseAcc = zeros(fftLen, 1);
    systemImpulse = zeros(fftLen - hop, 1);
    chirpAcc = zeros(fftLen, 1);
    chirpRec = zeros(fftLen - hop, 1);
    previousError = Inf;
    wndCorrectionWeightingPrev = wndCorrectionWeighting;
    for it = 1 : maxIter
        cnt = 1;
        systemImpulseAcc(:) = 0;
        systemImpulse(:) = 0;
        for i = 1 : (fftLen / hop) * 2
            if i <= (fftLen / hop) - 1
                stepSize = mod((fftLen - hop * i - fftLen / 2 - 1) - sht, fftLen);
                phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
                q_fft_frameinv = ltv_1st_2ndNoSq(phaseShifter, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, wndCorrectionWeighting, correctionWnd);
                q_fft_frameinv = q_fft_frameinv(:, 1);
                if reqSynthesisWnd
                    q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
                end
                q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
                yInvShifted = ifftshift(ifft(q_fft_frameinv));
                systemImpulseAcc(:) = systemImpulseAcc(:) + yInvShifted;
                myOut = systemImpulseAcc(1 : hop);
                systemImpulseAcc(:) = [systemImpulseAcc(hop + 1 : end); zeros(hop, 1)];
                systemImpulse(cnt : cnt + hop - 1) = myOut;
            end
            frame = bandlimitChirpPadded(frmIdx(i) : frmIdx(i) + fftLen - 1);
            X = fft(fftshift(frame));
            chirpSpectrum = ltv_1st_2ndNoSq(X(1 : halfLen), b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, wndCorrectionWeighting, correctionWnd);
            chirpSpectrum = chirpSpectrum(:, 1);
            % Virtually multiplying Hann window in time domain on frequency domain
            if reqSynthesisWnd
                chirpSpectrum = 0.25 * (2 * chirpSpectrum + [conj(chirpSpectrum(2)); chirpSpectrum(1 : end - 1)] + [chirpSpectrum(2 : end); conj(chirpSpectrum(end - 1))]);
            end
            % Reflect spectrum and conjugate
            chirpSpectrum(halfLen+1:fftLen,:) = conj(chirpSpectrum(halfLen-1:-1:2,:));
            chirpInvShifted = ifftshift(ifft(chirpSpectrum));
            % Overlap add
            chirpAcc(:) = chirpAcc(:) + chirpInvShifted;
            myOut = chirpAcc(1 : hop);
            chirpAcc(:) = [chirpAcc(hop + 1 : end); zeros(hop, 1)];
            chirpRec(cnt : cnt + hop - 1) = myOut;
            cnt = cnt + hop;
        end
        systemImpulse = circshift(systemImpulse, sht);
        systemImpulse(1 : sht) = 0;
        if hop ~= fftLen / 2
            truncatedSystemImpulse = [systemImpulse((fftLen - hop) - (fftLen/2) : end)];
        else
            truncatedSystemImpulse = [0; systemImpulse];
        end
        truncatedSystemImpulse = [truncatedSystemImpulse; truncatedSystemImpulse(halfLen - 1: -1 : 1)];
        chirpRec = circshift(chirpRec, -fftLen);
        errorCurr = sum(abs(chirpRec(1 : fftLen) - bandlimitChirp));
        if errorCurr >= previousError
            wndCorrectionWeighting = wndCorrectionWeightingPrev;
            break;
        else
            previousError = errorCurr;
            %         disp('iteration = ' + string(it) + ', error = ' + sprintf('%1.14f', errorCurr))
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
        SNR = 10*log10(sum(abs(target).^2)/sum(abs(target-synsq).^2)); %calculate the SNR
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
coeff.correctionWnd = correctionWnd;
coeff.wndCorrectionWeighting = wndCorrectionWeighting;
coeff.fftLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = hop;
coeff.fs = fs;
coeff.theoreticalWindowShape = theoreticalWindowShape;
end
function [f, g] = fitExp(corrF, cdf, raise)
f = mean( ( (corrF .^ raise) - cdf) .^ 2, 1); % Cost
J2 = -(2*corrF.^raise.*log(corrF) .* (cdf - corrF.^raise)); % Gradient
J2(isnan(J2)) = 0;
g = mean(J2, 1);
end
function x = chirp2(t,t1,f0,f1)
beta = (f1-f0)./t1;
x = cos(2*pi * ( 0.5* beta .* (t .* t) + f0 * t));
end
function spec = ltv_1st_2ndNoSq(dftSpec, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, wndCorrectionWeighting, correctinonWnd)
if ~isempty(wndCorrectionWeighting)
    q_fft_frameinv = dftSpec;
    q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    correctedTime = ifft(q_fft_frameinv);
    correctedTime2 = correctedTime .* correctinonWnd;
    getbackCorrectedToSpectrum2 = fft(correctedTime2);
    % Frequency weighted interpolation between original and windowed spectra
    dftSpec = dftSpec .* wndCorrectionWeighting + getbackCorrectedToSpectrum2(1 : halfLen) .* (1 - wndCorrectionWeighting);
end
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
spec = q_fft_frame(prepad + 1 : end - pospad + 1, :);
if ~isempty(corrF)
    spec = spec .* corrF;
end
spec(1, :) = real(spec(1, :));
spec(end, :) = real(spec(end, :));
end