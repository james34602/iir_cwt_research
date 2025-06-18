function [coeff, f] = IIR1(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
addpath('gradients/')
if nargin < 1
    fs = 48000;
    fftLen = 64;
    hop = 16;
    oct = 20;
    order = 2;
    reqSynthesisWnd = 1;
    sparCutOff = 5e-8;
    zp = 1;
end
outterWnd = 1;
paddedFFTLen = fftLen * zp;
eachSide = (paddedFFTLen - fftLen) / 2;
halfLen = getFFTHalfLen(fftLen);
paddedHalfLen = getFFTHalfLen(paddedFFTLen);
ovp = fftLen / hop;
% number of points of pre and post padding used to set initial conditions
prepad = fix(fftLen / 2);
pospad = fix(fftLen / 2);
% frequency bins grid (linear in this case) - pre and pos padding is added
% poles of the IIR LTV Q FFT transform for the parameters above
f = (0:1:fftLen/2)*fs/fftLen;
% number of points of pre and post padding used to set initial conditions
thetas = 0:(fftLen/2);
thetas(1) = eps;
%% Compute poles
sigmas = (thetas ./ fftLen) ./ oct / pi * fftLen;
if order == 2
    [b, a, c1, c2] = gauss_precompute(sigmas);
else
    a = 1 - exp(-1 ./ (0.3 / oct .* thetas.'));
    b = [];
    c1 = [];
    c2 = [];
end
%% Pole limiting
% analysisWnd = tukeywin_periodic(fftLen, 0.5);
analysisWnd = hann(fftLen, 'periodic') .^ (1 / reqSynthesisWnd);
synWnd1 = hann(fftLen, 'periodic') .^ reqSynthesisWnd;
if outterWnd
    synWnd2 = hann(paddedFFTLen, 'periodic') .^ reqSynthesisWnd;
    synWnd3 = synWnd2(eachSide + 1 : eachSide + fftLen);
else
    synWnd2 = [];
end
if zp == 1
    synWnd2(:) = synWnd1(:);
    synWnd3 = synWnd2(eachSide + 1 : eachSide + fftLen);
    synWnd1(:) = 1;
end
chopedWnd1 = analysisWnd(halfLen : end);
chopedWnd2 = synWnd1(halfLen : end);
halfWndLen = halfLen - 1;
digw2 = linspace(0, pi, halfLen);
digw = digw2(1 : halfWndLen);
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd1.';
h2 = h2 .* chopedWnd2.';
theoreticalWindowShape = [zeros(size(thetas, 2), 1), h2(:, (halfLen-1):-1:2), h2];
if outterWnd
    theoreticalWindowShape = theoreticalWindowShape .* synWnd3';
end
tol = min(1.5, max(1, (fftLen / hop) / 5));
hopsizeTol = min(fftLen, ceil(hop * 2 * tol));
if mod(hopsizeTol, 2) == 1
    hopsizeTol = hopsizeTol - 1;
end
smallestPossibleWnd = [zeros((fftLen - hopsizeTol) / 2, 1); hann(hopsizeTol, 'periodic'); zeros((fftLen - hopsizeTol) / 2, 1)];
wndDif = theoreticalWindowShape' - smallestPossibleWnd;
wndDifPwr = mean(abs(wndDif), 1);
[~, firstUndersampling1] = min(wndDifPwr);
firstUndersampling2 = find(any(wndDif < 0), 1, 'first');
firstUndersampling = ceil((firstUndersampling1 + firstUndersampling2) / 2);
firstUndersampling = max(firstUndersampling, fix(fftLen / hop * oct / 2));
firstUndersampling = [];
if ~isempty(firstUndersampling)
    thetas = 0:(fftLen/2);
    thetas(1) = eps;
    thetas(firstUndersampling : end) = thetas(firstUndersampling);
    thetas1 = thetas;
    thetas1(halfLen+1:fftLen) = conj(thetas1(halfLen-1:-1:2));
else
    thetas1 = thetas;
    thetas1(halfLen+1:fftLen) = conj(thetas1(halfLen-1:-1:2));
end
% Eliminate oscillation around corner
time = 0.026 * fftLen; % More elements in array less smoothing is needed
alpha = 1 / (1 + time);
[b2, a2] = butter(1, alpha);
thetas1 = filtfilt(b2, a2, thetas1);
mtheta = 10;
if min(thetas1) > mtheta
    flr = min(thetas1);
    cel = max(thetas1);
    thetas1 = thetas1 - flr;
    thetas1 = thetas1 ./ max(thetas1);
    thetas1 = thetas1 * (cel - mtheta) + mtheta;
end
thetas1(thetas1 < eps) = eps;
sigmas = (thetas1 ./ fftLen) ./ oct / pi * fftLen;
% sigmas(:) = max(sigmas);
if order == 2
    [b, a, c1, c2] = gauss_precompute(sigmas);
else
    a = 1 - exp(-1 ./ (0.3 / oct .* thetas1.'));
    b = [];
    c1 = [];
    c2 = [];
end
if any(abs(thetas1 - mean(thetas1)) < eps)
    flat = 1;
else
    flat = 0;
end
disp('First constant Q frequnecy that undersample ' + string(f(firstUndersampling)) + ' Hz')
%% Compute reassignment frequency domain sample shifter
phaseShifter1 = exp(1i * (2*pi*(0:halfLen-1)/fftLen).');
%% Obtain DFT filterbank frequency response
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd1.';
h2 = h2 .* chopedWnd2.';
theoreticalWindowShape = [zeros(size(thetas1, 2), 1), h2(:, (halfLen-1):-1:2), h2];
theoreticalWindowShape = [zeros(size(theoreticalWindowShape, 1), eachSide), theoreticalWindowShape, zeros(size(theoreticalWindowShape, 1), eachSide)];
%%
theoreticalWindowShape1 = theoreticalWindowShape(:, eachSide + 1 : eachSide + fftLen);
actualWindowShape3 = zeros(size(theoreticalWindowShape1));
dftMtx = dftmtx(fftLen);
Xa = fft(fftshift(analysisWnd));
% Xa = zeros(fftLen, 1);
% Xa(1) = fftLen/2;
% Xa(2) = fftLen/4;
% Xa(end) = fftLen/4;
% Xa3 = fft(dftMtx' .* fftshift(analysisWnd)).';
% plot(del2(b));
% hold on;
% plot(fftshift(del2(b)));hold off;
% axis tight
shtfft = -(2*rem(0:(halfLen-1), 2) - 1)';
for idx = 1 : size(theoreticalWindowShape1, 1)
    Xa2 = circshift(Xa, idx - 1);
    % Xa2_ = Xa3(:, idx);
    % Xa2 and Xa2_ is equal symbolically, but not numerically exact
    % disp(sum(abs(Xa2 - Xa2_)))
    yy = ltv_1st_2ndNoSq2(Xa2, b, a, c1, c2, prepad, pospad, [], halfLen);
    kk = fft(yy);
    b2 = fftshift(conj(dftMtx(:, idx)) .* kk);
    % plot(real(b2));
    % hold on
    % plot(imag(b2))
    % hold off
    % axis tight
    actualWindowShape3(idx, :) = b2;
end
actualWindowShape3 = (actualWindowShape3 .* synWnd2.') / fftLen;
windowedDFTMatrix2 = dftMtx .* actualWindowShape3;
LSsolution = ifft(windowedDFTMatrix2.');
LSsolution = real(LSsolution);
% kk = real(conj(dftMtx) .* fft(LSsolution)).';
Phi0 = LSsolution * dftMtx;
Phi0 = Phi0(1 : halfLen, :);

sigLen = 32768;
x = randn(sigLen, 1);
x = [zeros(fftLen * 0, 1); x(:); zeros(fftLen, 1)];
nframes = fix(ceil(length(x)/hop) - fftLen / hop / 2);
% zero padding at the end to complete the last frame
paddedLen = hop * nframes;
x = [x; zeros(length(x) - paddedLen + hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * hop;
y = zeros(nframes * hop, 1);
for m = 1 : nframes % index of window position
    Xa2 = LSsolution * fft(x(frmIdx(m) : frmIdx(m) + fftLen - 1)); % compute analysis filter outputs of a block of L samples
    Xa2 = Xa2(1 : halfLen);
    Xa = fft(fftshift(x(frmIdx(m) : frmIdx(m) + fftLen - 1) .* analysisWnd)); % compute analysis filter outputs of a block of L samples
    q_fft_frameinv = ltv_1st_2ndNoSq2(Xa, b, a, c1, c2, prepad, pospad, [], halfLen);
    q_fft_frameinv = 0.5 * (q_fft_frameinv + 0.5 * circshift(q_fft_frameinv, 1) + 0.5 * circshift(q_fft_frameinv, -1));
    q_fft_frameinv = shtfft .* q_fft_frameinv(1 : halfLen);
    plot(abs(Xa2-q_fft_frameinv));
    axis tight
    disp('')
end
%% Save coefficient
coeff.correctionMatrix = sparse(correctionMatrix);
coeff.prepad = prepad;
coeff.pospad = pospad;
coeff.b = b;
coeff.a = a;
coeff.c1 = c1;
coeff.c2 = c2;
coeff.phaseShifter1 = phaseShifter1;
coeff.correctionWnd = correctionWndHF;
coeff.fftLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = hop;
coeff.fs = fs;
coeff.reqSynthesisWnd = reqSynthesisWnd;
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2, prepad, pospad, corrF, halfLen)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end