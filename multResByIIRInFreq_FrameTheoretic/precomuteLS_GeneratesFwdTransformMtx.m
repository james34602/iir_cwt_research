function coeff = precomuteLS_GeneratesFwdTransformMtx(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 44100;
    fftLen = 2048;
    hop = 512;
    oct = 50;
    order = 2;
    reqSynthesisWnd = 1;
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
% analysisWnd(:)=1;
chopedWnd1 = analysisWnd(halfLen : end);
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
theoreticalWindowShape = [zeros(size(thetas, 2), 1), h2(:, (halfLen-1):-1:2), h2];
tol = min(1.5, max(1, (fftLen / hop) / 5));
hopsizeTol = min(fftLen, ceil(hop * 2 * tol));
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
    [b, a, c1, c2] = gauss_precompute(sigmas);
else
    a = 1 - exp(-1 ./ (0.3 / oct .* thetas1.'));
    b = [];
    c1 = [];
    c2 = [];
end
disp('First constant Q frequnecy that undersample ' + string(f(firstUndersampling)) + ' Hz')
origHalfLen = halfLen;
origFFTLen = fftLen;
actualWindowShape3 = zeros(fftLen, fftLen);
dftMtx = dftmtx(fftLen);
Xa = fft(fftshift(analysisWnd));
for idx = 1 : fftLen
    Xa2 = circshift(Xa, idx - 1);
    yy = ltv_1st_2ndNoSq2(Xa2, b, a, c1, c2);
    b2 = fftshift(conj(dftMtx(:, idx)) .* fft(yy));
    actualWindowShape3(idx, :) = b2;
end
windowedDFTMatrix1 = dftMtx .* actualWindowShape3;
LSsolution1 = ifft(windowedDFTMatrix1.') / fftLen;
Phi0 = LSsolution1 * dftMtx;
sig = randn(fftLen, 1);
trans1 = Phi0 * sig;
trans1 = trans1(1 : halfLen);
X = fft(fftshift(sig .* analysisWnd));
q_fft_frameinv = ltv_1st_2ndNoSq2(X, b, a, c1, c2);
shtfft2 = -(2*rem(0:(halfLen-1), 2) - 1)';
q_fft_frameinv2 = shtfft2 .* q_fft_frameinv(1 : origHalfLen);
%%
LSsolution1 = single(real(LSsolution1(1 : origHalfLen, :)));
q_fft_frameinv_ = LSsolution1 * fft(sig);
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
coeff.LSsolution1 = LSsolution1;
return;
outterPad = 0;
%%
[x, fs] = audioread('D:\DSP_research\time-frequency-analysis\iir_cqt_research\multResByIIRInFreq_FullyFreqDomainCorrMatrix\FanChirpF0gram_v1\pop1_long.wav');
% [x, fs] = loadSignal(7);
f = (0:1:fftLen/2)*fs/fftLen;
[specgram, dftspecgram, frmIdx, front, back] = ltv_spectrogram2(x, coeff, outterPad);
s_q = abs(specgram) .^ 2;
s_q = s_q ./ max(abs(s_q(:)));
dftspecgram = abs(dftspecgram) .^ 2;
dftspecgram = dftspecgram ./ max(abs(dftspecgram(:)));
R=renyi(s_q,frmIdx,f',3)
R2=renyi(dftspecgram,frmIdx,f',3)
%% Plot spectrogram
imagesc(mag2db(sqrt(s_q)))
caxis([-100, 0])
% ylim([100, 400])
% xlim([200, 400])
% colormap(jet);
black = abs(1-gray);
colormap(black);
set(gca,'YDir','normal');
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
    Sdft(i, :) = shtfft2 .* X(1 : coeff.origHalfLen);
    q_fft_frameinv = ltv_1st_2ndNoSq2(X, coeff.b, coeff.a, coeff.c1, coeff.c2);
    q_fft_frameinv2 = shtfft2 .* q_fft_frameinv(1 : coeff.origHalfLen);
    q_fft_frameinv2(1) = real(q_fft_frameinv2(1));
    q_fft_frameinv2(end) = real(q_fft_frameinv2(end));
    S(i, :) = q_fft_frameinv2;
end
S = permute(S, [2, 1, 3]);
Sdft = permute(Sdft, [2, 1, 3]);
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end
function [b, a, c1, c2] = gauss_precompute(sigma)
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
end