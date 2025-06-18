function [dftspecgram, specgram, specgram4th, f, frmIdx, coeff] = HighOrderLS_Spectrogram(x, fftLen, hop, fs, oct, order, reqSynthesisWnd, STFTanalysisWnd)
rng(1)
% addpath('D:\DSP_research\time-frequency-analysis\iir_cqt_research\SVF_test')
if nargin < 1
    fs = 48000;
    fftLen = 4096;
    hop = 128;
    oct = 50;
    order = 2;
    reqSynthesisWnd = 2;
end
halfLen = getFFTHalfLen(fftLen);
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
sigmas = thetas / (oct * pi);
if order == 2
    [b, a, c1, c2] = gauss_precompute(sigmas);
else
    a = 1 - exp(-1 ./ (0.3 / oct .* thetas.'));
    b = [];
    c1 = [];
    c2 = [];
end
%% Pole limiting
analysisWnd = hann(fftLen, 'periodic') .^ (1 / reqSynthesisWnd);
synWnd1 = hann(fftLen, 'periodic') .^ reqSynthesisWnd;
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
tol = min(1.5, max(1, (fftLen / hop) / 5));
hopsizeTol = min(fftLen, ceil(hop * tol));
if mod(hopsizeTol, 2) == 1
    hopsizeTol = hopsizeTol - 1;
end
smallestPossibleWnd = [zeros((fftLen - hopsizeTol) / 2, 1); hann(hopsizeTol, 'periodic'); zeros((fftLen - hopsizeTol) / 2, 1)];
wndDif = theoreticalWindowShape' - smallestPossibleWnd;
wndDifPwr = sum(abs(wndDif), 1);
[~, firstUndersampling1] = min(wndDifPwr);
firstUndersampling2 = find(any(wndDif < 0), 1, 'first');
firstUndersampling = ceil((firstUndersampling1 + firstUndersampling2) / 2);
% firstUndersampling = max(firstUndersampling, fix(fftLen / hop * oct / 2));
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
    [b, a, c1, c2, divareaSum2nd] = gauss_precompute2(sigmas, length(sigmas));
    coeff.divareaSum2nd = divareaSum2nd;
    [~, ~, ~, Bmfwd, Bmbwd, Am, ~, ~, ~] = InitDeriche(sigmas(1), length(sigmas));
    Bmfwd = zeros([size(Bmfwd), length(sigmas)]);
    Bmbwd = zeros([size(Bmbwd), length(sigmas)]);
    Am = zeros([size(Am), length(sigmas)]);
    divareaSum = zeros(length(sigmas), 1);
    for i = 1 : length(sigmas)
        [~, ~, ~, Bmfwd(:, :, i), Bmbwd(:, :, i), Am(:, :, i), ~, ~, divareaSum(i)] = InitDeriche(sigmas(i), length(sigmas));
    end
    coeff.Bmfwd = Bmfwd;
    coeff.Bmbwd = Bmbwd;
    coeff.Am = Am;
    coeff.divareaSum = divareaSum;
else
    a = 1 - exp(-1 ./ (0.3 / oct .* thetas1.'));
    b = [];
    c1 = [];
    c2 = [];
end
disp('First constant Q frequnecy that undersample ' + string(f(firstUndersampling)) + ' Hz')
%% Compute reassignment frequency domain sample shifter
% phaseShifter1 = exp(1i * (2*pi*(0:halfLen-1)/fftLen).');
origHalfLen = halfLen;
origFFTLen = fftLen;
%%
coeff.origHalfLen = origHalfLen;
coeff.origFFTLen = origFFTLen;
coeff.hop = hop;
coeff.prepad = prepad;
coeff.pospad = pospad;
coeff.b = b;
coeff.a = a;
coeff.c1 = c1;
coeff.c2 = c2;
coeff.analysisWnd = analysisWnd;
coeff.STFTanalysisWnd = STFTanalysisWnd;
outterPad = 0;
coeff.corrF = ones(halfLen, 1);
coeff.f = f;
%%
f = (0:1:fftLen/2)*fs/fftLen;
[specgram4th, ~, ~, ~] = ltv_spectrogram1(x, coeff, outterPad);
[specgram, dftspecgram, frmIdx, front, back] = ltv_spectrogram2(x, coeff, outterPad);
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function [S, Sdft, frmIdx, front, back] = ltv_spectrogram2(x, coeff, outterPad)
x = x(:);
xLen = length(x);
x = [zeros(outterPad + coeff.origFFTLen - coeff.hop, 1); x; zeros(outterPad + coeff.origFFTLen / 2, 1)];
front = outterPad + coeff.origFFTLen - coeff.hop;
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.origFFTLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
tmp = length(x) - paddedLen + coeff.hop * 2;
x = [x; zeros(tmp, 1)];
back = length(x) - front - xLen;
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
halfLen = coeff.origFFTLen/2+1;
shtfft2 = -(2*rem(0:(halfLen-1), 2) - 1)';
S = zeros(nframes, halfLen, 'like', 1i);
Sdft = zeros(nframes, halfLen, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.origFFTLen - 1);
    X = fft(fftshift(frame .* coeff.analysisWnd));
    X2 = fft(fftshift(frame .* coeff.STFTanalysisWnd));
    Sdft(i, :) = shtfft2 .* X2(1 : coeff.origHalfLen);
    q_fft_frameinv = ltv_1st_2ndNoSq2(X, coeff.b, coeff.a, coeff.c1, coeff.c2);
    q_fft_frameinv2 = shtfft2 .* q_fft_frameinv(1 : coeff.origHalfLen);
    q_fft_frameinv2(1) = real(q_fft_frameinv2(1));
    q_fft_frameinv2(end) = real(q_fft_frameinv2(end));
    S(i, :) = q_fft_frameinv2;
end
S = permute(S, [2, 1, 3]);
Sdft = permute(Sdft, [2, 1, 3]);
frmIdx = frmIdx - coeff.origHalfLen + coeff.hop;
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end
function [S, frmIdx, front, back] = ltv_spectrogram1(x, coeff, outterPad)
x = x(:);
xLen = length(x);
x = [zeros(outterPad + coeff.origFFTLen - coeff.hop, 1); x; zeros(outterPad + coeff.origFFTLen / 2, 1)];
front = outterPad + coeff.origFFTLen - coeff.hop;
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.origFFTLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
tmp = length(x) - paddedLen + coeff.hop * 2;
x = [x; zeros(tmp, 1)];
back = length(x) - front - xLen;
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
halfLen = coeff.origFFTLen/2+1;
shtfft2 = -(2*rem(0:(halfLen-1), 2) - 1)';
S = zeros(nframes, halfLen, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.origFFTLen - 1);
    X = fft(fftshift(frame .* coeff.analysisWnd));
    q_fft_frameinv = ltv_1st_2ndNoSq1(X, coeff.Bmfwd, coeff.Bmbwd, coeff.Am, coeff.divareaSum);
    q_fft_frameinv2 = shtfft2 .* q_fft_frameinv(1 : coeff.origHalfLen);
    q_fft_frameinv2(1) = real(q_fft_frameinv2(1));
    q_fft_frameinv2(end) = real(q_fft_frameinv2(end));
    S(i, :) = q_fft_frameinv2;
end
S = permute(S, [2, 1, 3]);
frmIdx = frmIdx - coeff.origHalfLen + coeff.hop;
end
function spec = ltv_1st_2ndNoSq1(dftSpec, Bmfwd, Bmbwd, Am, divareaSum)
%% Gaussian windowing
spec = smpparsvffiltbidirectional2(Bmfwd, Bmbwd, Am, dftSpec) .* divareaSum;
end
function y = smpparsvffiltbidirectional2(b1, b2, a, x)
c1 = a(2, :, :) + 2;
c2 = (1 + a(2, :, :) + a(3, :, :)) ./ c1;
d0fwd = b1(1, :, :);
d1fwd = (2 * b1(1, :, :) + b1(2, :, :)) ./ c1;
d2fwd = (b1(1, :, :) + b1(2, :, :)) ./ (c1 .* c2);
d0bwd = b2(1, :, :);
d1bwd = (2 * b2(1, :, :) + b2(2, :, :)) ./ c1;
d2bwd = (b2(1, :, :) + b2(2, :, :)) ./ (c1 .* c2);
z1_A = [0, 0];
z2_A = [0, 0];
y1 = zeros(size(x));
for i = 1 : length(x)
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0fwd(1, 1, i) * Y + d1fwd(1, 1, i) * z1_A(1) + d2fwd(1, 1, i) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0fwd(1, 2, i) * Y2 + d1fwd(1, 2, i) * z1_A(2) + d2fwd(1, 2, i) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i) * Y2;
    y1(i) = Yi1 + Yi2;
end
z1_A = [0, 0];
z2_A = [0, 0];
y2 = zeros(size(x, 1), 1);
for i = length(x) : -1 : 2
    Xi = x(i);
    Y = Xi - z1_A(1) - z2_A(1);
    Yi1 = d0bwd(1, 1, i - 1) * Y + d1bwd(1, 1, i - 1) * z1_A(1) + d2bwd(1, 1, i - 1) * z2_A(1);
    z2_A(1) = z2_A(1) + c2(1, 1, i - 1) * z1_A(1);
    z1_A(1) = z1_A(1) + c1(1, 1, i - 1) * Y;
    Y2 = Xi - z1_A(2) - z2_A(2);
    Yi2 = d0bwd(1, 2, i - 1) * Y2 + d1bwd(1, 2, i - 1) * z1_A(2) + d2bwd(1, 2, i - 1) * z2_A(2);
    z2_A(2) = z2_A(2) + c2(1, 2, i - 1) * z1_A(2);
    z1_A(2) = z1_A(2) + c1(1, 2, i - 1) * Y2;
    y2(i - 1) = Yi1 + Yi2;
end
y = y1 + y2;
end
function [b, a, c1, c2, divareaSum] = gauss_precompute2(sigma, intendedArrayLen)
sigma = sigma(:);
fnc1 = @(x) ((606 * x.^2) / 1087 - (3009 * x) / 5513 + 712 / 5411) ./ (x.^2 - (496 * x) / 541 + 719 / 1034);
mp = exp(-(137 / 100) ./ sigma);
a1 = -2 * cos(fnc1(sigma) ./ sigma) .* mp;
%% Transfer function
a = [ones(size(sigma, 1), 1), a1, mp .* mp];
b = sum(a, 2);
%% SVF realization
c1 = 2 + a1;
c2 = b ./ c1;
x = zeros(intendedArrayLen, 1);
x(getFFTHalfLen(intendedArrayLen)) = 1;
divareaSum = zeros(size(b));
for idx = 1 : length(sigma)
    yp1 = filter(b(idx), a(idx, :), x);
    ym1 = flipud(filter(b(idx), a(idx, :), flipud(yp1)));
    ym1 = ym1 / max(ym1);
    divareaSum(idx) = 1 / sum(ym1);
end
end