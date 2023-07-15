function coeff = ltv_precomute2(fftLen, fs, oct, order)
addpath('minFunc_2012/autoDif')
addpath('minFunc_2012/minFunc')
addpath('minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 512;
    oct = 4;
    order = 2;
end
zp = 1;
halfLen = getFFTHalfLen(fftLen);
% number of points of pre and post padding used to set initial conditions
prepad = 10;
pospad = 42;
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
%plot(ifftshift(ifft(q_fft_frame)) .* 1)
%disp(sigmas(prepad+sel))
%% Compute reassignment frequency domain sample shifter
phaseShifter1 = exp(1i * (2*pi*(0:halfLen-1)/fftLen).');
phaseShifter2 = exp(-1i * (2*pi/fftLen).');
phaseShifter3 = exp(1i * (2*pi/fftLen).');
%% Obtain filterbank overall response
tmp = zeros(fftLen, fftLen + 1);
tmp2 = zeros(halfLen, fftLen + 1);
corrF = [];
for i = 0 : fftLen
    stepSize = mod(fftLen - i - fftLen / 2 - 1, fftLen);
    phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
    q_fft_frameinv = ltv_1st_2ndNoSq2(phaseShifter, b, a, c1, c2, prepad, pospad, corrF, halfLen);
    %% Inverse transform
    tmp2(:, i + 1) = q_fft_frameinv;
    % Reflect spectrum and conjugate
    q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    yInvShifted = ifftshift(ifft(q_fft_frameinv)).';
    tmp(:, i + 1) = yInvShifted;
end
systemImpulse = overlapAdd(tmp, 1);
systemImpulse = systemImpulse((halfLen : halfLen + fftLen - 1) - 1);
% systemImpulse2 = sum(tmp2(:, :), 1) + conj( sum(tmp2(2 : fftLen / 2, :), 1) );
truncatedSystemImpulse = fft(systemImpulse);
truncatedSystemImpulse = abs(truncatedSystemImpulse(1 : halfLen));
corrF = 1 ./ truncatedSystemImpulse;
%% Pole limiting
cornerSmoothing = 1;
gaussSigmaLimiting = 1;
firstUndersampling = floor(halfLen * 0.8);
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
%% Compute correction for synchrosqueezing transform
corrS = ones(halfLen, 1);
max1 = -Inf;
firstBestSol = zeros(halfLen, 1);
it = 0;
target = zeros(fftLen, 1);
target(fftLen / 2) = 1;
kDeltaSpectrum = zeros(halfLen, fftLen);
for i=1:fftLen
    stepSize = mod(fftLen - i - fftLen / 2, fftLen);
    phaseShifter = exp(-1i * stepSize * (2*pi*(0:halfLen-1)/fftLen)');
    S2 = ltv_1st_2nd(phaseShifter, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, 1);
    kDeltaSpectrum(:, i) = real(S2(:, 1));
end
while (1)
    S2 = kDeltaSpectrum .* corrS;
    synsq = sum(S2, 1) + conj( sum(S2(2 : fftLen / 2, :), 1) );
    synsq = synsq(:);
    synsqSpec = fft(synsq);
    synsqSpec = abs(synsqSpec(1 : halfLen));
    SNR = 10*log10(1 ./ sum(abs(target-synsq).^2)); % Calculate the SNR
    it = it + 1;
    if it > 100
        break;
    end
    if SNR > max1
        max1 = SNR;
        firstBestSol(:) = corrS;
    end
    corrS = corrS .* 1 ./ synsqSpec;
end
S2 = kDeltaSpectrum .* firstBestSol;
synsq = sum(S2(:, :), 1) + conj( sum(S2(2 : fftLen / 2, :), 1) );
synsq = real(synsq(:));
%% Optimization that doesn't works end
iv = ifft(1 ./ fft(synsq));
%% Save coefficient
coeff.ir = iv;
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
coeff.fftLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = 1;
coeff.fs = fs;
coeff.f = (0 : (halfLen - 1)) / fftLen * fs;
coeff.k = coeff.f / fs * fftLen;
coeff.W = exp(coeff.k * (2 * pi) / fftLen * 1i);
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2, prepad, pospad, corrF, halfLen)
%% Hann in frequency domain with input being shifted by fftshift
specHann = 2 * dftSpec + [conj(dftSpec(2)); dftSpec(1 : end - 1)] + [dftSpec(2 : end); conj(dftSpec(end - 1))];
% specHann = 2 * specHann + [conj(specHann(2)); specHann(1 : end - 1)] + [specHann(2 : end); conj(specHann(end - 1))];
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
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end