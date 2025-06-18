function [coeff, f] = precomuteLS_LTI_test(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 512;
    hop = 32;
    oct = 20;
    order = 2;
    reqSynthesisWnd = 1;
    sparCutOff = 5e-8;
    zp = 1;
end
ovp = fftLen / hop;
paddedFFTLen = fftLen * zp;
eachSide = (paddedFFTLen - fftLen) / 2;
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
% Xa3 = fft(diag(fftshift(analysisWnd)) * dftMtx'); % Matrix multiplication representation
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
    yy = ltv_1st_2ndNoSq2(Xa2, b, a, c1, c2);
    b2 = fftshift(conj(dftMtx(:, idx)) .* fft(yy));
    actualWindowShape3(idx, :) = b2;
end
actualWindowShape3 = (actualWindowShape3 .* synWnd1.') / fftLen;
%% Plot result after overlap-add
clear cpxRes dftFilterbank2 theoreticalWindowShape2 theoreticalWindowShape
windowedDFTMatrix1 = dftMtx .* actualWindowShape3;
LSsolution1 = ifft(windowedDFTMatrix1.');
Phi0 = LSsolution1 * dftMtx;
Phi0 = Phi0(1 : halfLen, :);
Phi0(1, :) = real(Phi0(1, :));
Phi0(end, :) = real(Phi0(end, :));
Phi0(halfLen+1:fftLen, :) = conj(Phi0(halfLen-1:-1:2, :));
origHalfLen = halfLen;
origFFTLen = fftLen;
[halfLen, fftLen] = size(Phi0);
% build frame matrix Phi of Xa = Phi*x
% Phi is used to compute pseudo inverse (dual frame matrix)
ckmin = 3;
% ckmin = 2;
ck = max(ckmin,floor(fftLen/halfLen)); % determines size of matrix Phi (and reconstruction error)
Nx = ck*fftLen;
Nm = ceil(Nx/hop);
k = ck*ceil(fftLen/hop);
theoreticalCentre = ceil(k*halfLen/(2*halfLen)) * hop + origHalfLen - fftLen;
offset = getFFTHalfLen(ck*fftLen-fftLen) - theoreticalCentre;
if hop*(Nm-1) + fftLen < ck*fftLen-offset
    disp([hop*(Nm-1) + fftLen, ck*fftLen-offset])
    error('Synthesis frame matrix does''t exist');
end
tic
tPhi0 = Phi0';
blk = real(tPhi0*Phi0);
w = size(blk,2);
h = size(blk,1);
N = Nm;
vshft = hop;
hshft = hop;
nRows = vshft*(N-1) + h;
nCols = hshft*(N-1) + w;
C = sparse(nRows,nCols);
[jdx0,idx0] = meshgrid(1:w,1:h);
for l = 1:N
    idx = idx0 + (l-1)*vshft;
    jdx = jdx0 + (l-1)*hshft;
    if l > 2
        C = C + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    else
        C = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    end
end
toc
%%
epsi = 1e-14;
C2 = C(fftLen-offset+1:ck*fftLen-offset, (fftLen-offset+1:ck*fftLen-offset)) + epsi*speye(ck*fftLen-fftLen);
% imagesc(log(abs(C2)))
padUp = ceil(k*halfLen/(2*halfLen)) * hop + offset;
padDown = hop*(Nm-1) - padUp;
tPhiSel = [zeros(padUp, halfLen); tPhi0; zeros(padDown, halfLen)];
cuttedPos = tPhiSel(fftLen+1:ck*fftLen, :);
Phi_s2 = C2 \ cuttedPos; % pseudo inverse
clear C C2 tPhiSel
Ls = size(Phi_s2,1);
newHalfLen = getFFTHalfLen(size(Phi_s2, 1));
% dmt = dftmtx(size(Phi_s2, 1))';
% sol2 = real(dmt\Phi_s2 * size(Phi_s2, 1));
sol2 = real(fft(Phi_s2, [], 1));
halfLen2 = size(Phi_s2, 1) / 2 + 1;
sol3 = sol2(1 : halfLen2, :);
sparCutOff = 5e-15;
gWeights = abs(sol3) > sparCutOff;
sparsity = 1 - sum(gWeights(:)) / numel(sol3);
sol3=sparse(sol3.*gWeights);
sigLen = 32768;
x = [zeros(fftLen, 1); randn(sigLen, 1); zeros(fftLen, 1)];
sht = 10;
x2 = circshift(x, sht);
coeff.origHalfLen = origHalfLen;
coeff.origFFTLen = origFFTLen;
coeff.Ls = Ls;
coeff.hop = hop;
coeff.prepad = prepad;
coeff.pospad = pospad;
coeff.b = b;
coeff.a = a;
coeff.c1 = c1;
coeff.c2 = c2;
coeff.analysisWnd = analysisWnd;
outterPad = 0;
%%
specgram1 = ltv_spectrogram2(x, coeff, outterPad);
specgram1(origHalfLen+1:origFFTLen, :) = conj(specgram1(origHalfLen-1:-1:2, :));
E1 = norm(abs(specgram1),'fro').^2;
sg.SRe = real(specgram1);
sg.SIm = imag(specgram1);
sg.target = x;
identity = full(sol3(:));
func1 = @(x) corrRealMtxFastExtendedPagedHC2(x, origHalfLen, origFFTLen, newHalfLen, Ls, hop, sg, [], outterPad);
specgram2 = ltv_spectrogram2(x2, coeff, outterPad);
specgram2(origHalfLen+1:origFFTLen, :) = conj(specgram2(origHalfLen-1:-1:2, :));
E2 = norm(abs(specgram2),'fro').^2;
disp(E2 - E1)
sg.SRe = real(specgram2);
sg.SIm = imag(specgram2);
sg.target = x2;
func2 = @(x) corrRealMtxFastExtendedPagedHC2(x, origHalfLen, origFFTLen, newHalfLen, Ls, hop, sg, [], outterPad);
[fval, grad, y1] = func1(identity);
[fval, grad, y2] = func2(identity);
y2Rec = circshift(y2, -sht);
dif = y1 - y2Rec;
difMAE = sum(abs(dif(fftLen + 1 : end - fftLen)))
end
function [f, wt_grad, y] = corrRealMtxFastExtendedPagedHC2(weights, origHalfLen, origFFTLen, newHalfLen, Ls, hop, sigs, gWeights, outterPad)
extraLatency = newHalfLen - origHalfLen;
x = reshape(weights, newHalfLen, origFFTLen);
q_fft_frameinv2Re = pagemtimes(x, sigs.SRe);
q_fft_frameinv2Im = pagemtimes(x, sigs.SIm);
q_fft_frameinv2Im(1, :, :) = 0;
q_fft_frameinv2Im(end, :, :) = 0;
if mod(Ls,2) == 0
    fg = 0;
    q_fft_frameinv2Re(newHalfLen+1:Ls, :, :) = q_fft_frameinv2Re(newHalfLen-1:-1:2, :, :);
    q_fft_frameinv2Im(newHalfLen+1:Ls, :, :) = -q_fft_frameinv2Im(newHalfLen-1:-1:2, :, :);
else
    fg = 1;
    q_fft_frameinv2Re(newHalfLen+1:Ls, :, :) = q_fft_frameinv2Re(newHalfLen:-1:2, :, :);
    q_fft_frameinv2Im(newHalfLen+1:Ls, :, :) = -q_fft_frameinv2Im(newHalfLen:-1:2, :, :);
end
mtx = ifft(q_fft_frameinv2Re + q_fft_frameinv2Im * 1j);
y3 = overlapAdd3D(mtx, hop);
y = y3(outterPad + origFFTLen - hop + extraLatency + 1 : outterPad + extraLatency + size(sigs.target, 1) - hop + origFFTLen, :);
dif = y - sigs.target;
avgG = 1 / (size(sigs.target, 1) * size(sigs.target, 2));
f = sum(dif(:) .^ 2) * avgG;
se_grad = ones(size(sigs.target)) * avgG .* dif * 2;
y_grad = zeros(size(y3));
y_grad(outterPad + origFFTLen - hop + extraLatency + 1 : outterPad + extraLatency + size(sigs.target, 1) - hop + origFFTLen, :) = se_grad;
idealBufferOutLen = ceil((size(y_grad, 1)-(Ls - hop))/(Ls-(Ls - hop)));
requiredBufferOutLen = size(sigs.SRe, 2);
cutLen = ceil(size(y_grad, 1) / hop) - requiredBufferOutLen;
rec = zeros(Ls, ceil(size(y_grad, 1) / hop), size(y_grad, 2));
for idx = 1 : size(y_grad, 2)
    rec(:, :, idx) = buffer(y_grad(:, idx), Ls, Ls - hop);
end
% tmp = ifft(rec(:, Ls / hop : end, :));
tmp = ifft(rec(:, cutLen + 1 : end, :));
if fg == 0
    tmp = tmp(1 : newHalfLen, :, :);
    tmp(2 : end - 1, :, :) = tmp(2 : end - 1, :, :) * 2;
else
    tmp = tmp(1 : newHalfLen, :, :);
    tmp(2 : end, :, :) = tmp(2 : end, :, :) * 2;
end
xRe_grad = real(tmp);
xIm_grad = -imag(tmp);
xIm_grad(1, :, :) = 0;
xIm_grad(end, :, :) = 0;
wt_grad = sum(pagemtimes(xIm_grad, 'none', sigs.SIm, 'transpose'), 3) + sum(pagemtimes(xRe_grad, 'none', sigs.SRe, 'transpose'), 3);
wt_grad = wt_grad(:);
if ~isempty(gWeights)
    wt_grad = wt_grad .* gWeights;
end
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function S = ltv_spectrogram2(x, coeff, outterPad)
x = x(:);
x = [zeros(outterPad + coeff.origFFTLen - coeff.hop, 1); x; zeros(outterPad + coeff.origFFTLen / 2, 1)];
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.origFFTLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
x = [x; zeros(length(x) - paddedLen + coeff.hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
halfLen = coeff.origFFTLen/2+1;
shtfft2 = -(2*rem(0:(halfLen-1), 2) - 1)';
S = zeros(nframes, halfLen, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.origFFTLen - 1);
    X = fft(fftshift(frame .* coeff.analysisWnd));
    q_fft_frameinv = ltv_1st_2ndNoSq2(X, coeff.b, coeff.a, coeff.c1, coeff.c2);
    lastElement = q_fft_frameinv(end);
    q_fft_frameinv = q_fft_frameinv(1 : coeff.origHalfLen + 1);
    q_fft_frameinv = 0.5 * (q_fft_frameinv + 0.5 * [lastElement; q_fft_frameinv(1 : end - 1)] + 0.5 * [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
    q_fft_frameinv2 = shtfft2 .* q_fft_frameinv(1 : coeff.origHalfLen);
    q_fft_frameinv2(1) = real(q_fft_frameinv2(1));
    q_fft_frameinv2(end) = real(q_fft_frameinv2(end));
    S(i, :) = q_fft_frameinv2;
end
S = permute(S, [2, 1, 3]);
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end
function [L, k] = cholesky(A, effectiveLen, fftLen, hop, Nm)
[~, n] = size(A);
L = zeros(size(A));
outLen = fftLen - mod((1 : n) - 1, hop);
for k = 1 : n
    nblk = k / hop;
    v1 = L(k : n, 1 : (k - 1));
    v2 = L(k, 1 : (k - 1))';
    L(k : n, k) = A(k : n, k) - v1 * v2;
    L(k, k) = sqrt(L(k, k));
    if (k < n)
        L((k + 1) : n, k) = L((k + 1) : n, k) / L(k, k);
    end
end
end