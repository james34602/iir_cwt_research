function [coeff, f] = precomuteLS_SparseGrad(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 64;
    hop = 10;
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
dbg2 = 1;
if dbg2
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
    % tic
    % Phi_ = sparse(halfLen*Nm,(Nm-1)*hop+fftLen);
    % for m = 0:Nm-1
    %     Phi_(m*halfLen+1:m*halfLen+halfLen,m*hop+1:m*hop+fftLen) = Phi0;
    % end
    % tPhi_ = Phi_';
    % C_ = real(tPhi_*Phi_);
    % toc
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
    N = Nm + 1;
    vshft = hop;
    hshft = hop;
    nRows = vshft*(N-1) + h;
    nCols = hshft*(N-1) + w;
    C3 = sparse(nRows,nCols);
    [jdx0,idx0] = meshgrid(1:w,1:h);
    for l = 1:N
        idx = idx0 + (l-1)*vshft;
        jdx = jdx0 + (l-1)*hshft;
        if l > 2
            C3 = C3 + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        else
            C3 = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        end
    end
    N = ceil((fftLen-offset)/hop) + ceil(fftLen/hop);
    vshft = hop;
    hshft = hop;
    nRows = vshft*(N-1) + h;
    nCols = hshft*(N-1) + w;
    pruned = zeros(nRows,nCols);
    [jdx0,idx0] = meshgrid(1:w,1:h);
    for l = 1:N
        idx = idx0 + (l-1)*vshft;
        jdx = jdx0 + (l-1)*hshft;
        if l > 2
            pruned = pruned + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        else
            pruned = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        end
    end
    epsi = 1e-14;
    pruned = full(pruned(fftLen-offset +1:fftLen-offset + fftLen+hop, fftLen-offset +1:fftLen-offset + fftLen+hop));
    pruned = pruned + epsi*speye(size(pruned));
    % effectiveLen = zeros(fftLen, 1);
    % prunedShifted = zeros(size(pruned));
    % for i = 1 : fftLen
    %     effectiveLen(i) = fftLen - i + 1;
    %     prunedShifted(:, i) = circshift(pruned(:, i), -(i - 1));
    %     prunedShifted(fftLen - i + 2 : end, i) = 0;
    % end
    % C2 = C(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen)) + epsi*speye(ck*fftLen-fftLen);
    % mtxeig = eig(full(C2));
    % frameBoundA_ = min(mtxeig)
    % frameBoundB_ = max(mtxeig)
    % pruned2 = C2(hop + 1 : hop + fftLen, hop + 1 : hop + fftLen);
    % R = cholesky(C2);
    % tt = tril(C2);
    % R2 = cholesky(tt, effectiveLen, fftLen, hop);
    toc
    fprintf(1,'computing inverse analysis frame matrix...\n');
    %%
    C=full(C);
    C2 = C(fftLen-offset+1:ck*fftLen-offset, (fftLen-offset+1:ck*fftLen-offset)) + epsi*speye(ck*fftLen-fftLen);
    C2_ = C3(fftLen-offset+1:ck*fftLen-offset, (fftLen-offset+1:ck*fftLen-offset)) + epsi*speye(ck*fftLen-fftLen);
    % C2 = C2_;
    % figure(1)
    % imagesc(log(abs(C2)))
    % figure(2)
    % imagesc(log(abs(pruned)))
    C4 = C(hop*floor((fftLen-offset)/hop) + 1: hop*ceil((ck*fftLen-offset)/hop), hop*floor((fftLen-offset)/hop) + 1: hop*ceil((ck*fftLen-offset)/hop));
    C4 = C4 + epsi*speye(size(C4));
    % imagesc(log(abs(C3)))
    pd = fftLen-offset - (hop*floor((fftLen-offset)/hop));
    % imagesc(log(abs(C2)))
    padUp = ceil(k*halfLen/(2*halfLen)) * hop + offset;
    padDown = hop*(Nm-1) - padUp;
    tPhiSel = [zeros(padUp, halfLen); tPhi0; zeros(padDown, halfLen)];
    % padUp = ic / halfLen * hop - fftLen;
    % if padUp < 1
    %     cut1 = (ck-1)*fftLen - fftLen - padUp;
    % else
    %     cut1 = 0;
    % end
    % cut2 = (ck-1)*fftLen - fftLen - padUp;
    % cc = [zeros(padUp, halfLen); tPhi0(cut1 + 1 : end, :); zeros(cut2, halfLen)];
    % disp(sum(sum(abs(tPhiSel(fftLen+1:ck*fftLen, :)-cc))))
    % cuttedPos = tPhiSel(fftLen+1:ck*fftLen, :);
    cuttedPos = tPhiSel(fftLen+1:ck*fftLen, :);
    Phi_s2 = C2 \ cuttedPos; % pseudo inverse
    Phi_s2_ = C2_ \ cuttedPos; % pseudo inverse
    trim = ((hop*ceil((ck*fftLen-offset)/hop)) - (hop*floor((fftLen-offset)/hop))) - (((ck*fftLen)-fftLen) + pd);
    cuttedPos2 = [zeros(pd, halfLen); cuttedPos; zeros(trim, halfLen)];
    %% Cholesky
    effectiveLen = zeros(fftLen, 1);
    tt2 = full(tril(C4));
    R2 = cholesky(tt2, effectiveLen, fftLen, hop, Nm);
    Phi_s4 = R2' \ (R2 \ cuttedPos2); % pseudo inverse
    Phi_s3 = Phi_s4(pd + 1 : pd + ((ck*fftLen)-fftLen), :);
    clear C C2 tPhiSel
else
    load('matlab2.mat');
end
% limit Phi_s2 according to chosen reconstruction error
tol = eps*15;
tol = tol*max(abs((Phi_s2(:))));
synMtxPwr = abs(Phi_s2) > tol;
imin = size(Phi_s2, 1);
imax = 1;
for idx = 1 : size(Phi_s2, 2)
    ii1 = find(synMtxPwr(:, idx), 1, 'first');
    ii2 = find(synMtxPwr(:, idx), 1, 'last');
    if ii1 < imin
        imin = ii1;
    end
    if ii2 > imax
        imax = ii2;
    end
end
imax = size(Phi_s2,1) - min(imin, size(Phi_s2,1) - imax);
if mod(imax-imin, 2) == 0
    imax = imax + 1;
end
% imin=imin-2;
% imax=imax-2;
Phi_s3 = Phi_s2(imin:imax,:); % Reduce latency
[~, mostMaxIdx] = max(Phi_s3(:, 1), [], 1);
newHalfLen = getFFTHalfLen(size(Phi_s3, 1));
offset = newHalfLen - mostMaxIdx;
if imin - offset > 0 && imax - offset < size(Phi_s2, 1)
    if mostMaxIdx ~= newHalfLen
        Phi_s3 = Phi_s2(imin - offset : imax - offset, :);
    end
    Phi_s2 = Phi_s3;
end
clear Phi_s3
Ls = size(Phi_s2,1);
newHalfLen = getFFTHalfLen(size(Phi_s2, 1));
% dmt = dftmtx(size(Phi_s2, 1))';
% sol2 = real(dmt\Phi_s2 * size(Phi_s2, 1));
sol2 = real(fft(Phi_s2, [], 1));
halfLen2 = size(Phi_s2, 1) / 2 + 1;
sol3 = sol2(1 : halfLen2, :);
sparCutOff = 5e-15;
gWeights = abs(sol3) > sparCutOff;
sparsity = 1 - sum(gWeights(:)) / numel(sol3)
sol3=sparse(sol3.*gWeights);
disp('')
% R = cholesky(C2);
% figure(1)
% imagesc(log10(abs(R) + eps))
% R2 = cholesky(pruned);
% figure(2)
% imagesc(log10(abs(R2 + eps)))
% R3 = R2(1 : fftLen + hop * (N - 2), 1 : fftLen + hop * (N - 2));
% figure(3)
% imagesc(log10(abs(R3 + eps)))
% R2 = chol(C2);
% Phi_s2 = R2 \ (R2' \ tPhi);
% res = backwardSubstitution(R', forwardSubstitution(R, full(tPhi(:, 2))));
% since Phi is complex
% fresp = fft(Phi_s2, [], 1);
% aa = fft(conj(dftMtx), [], 1);
% nDftMtx = dftmtx(Ls);
% nv = conj(nDftMtx) .* fresp;
fprintf(1,'size of synthesis frame matrix = %dx%d\n',Ls,halfLen);
sigLen = 32768;
x = randn(sigLen, 1);
%%
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
specgram = ltv_spectrogram2(x, coeff, outterPad);
specgram(origHalfLen+1:origFFTLen, :) = conj(specgram(origHalfLen-1:-1:2, :));
sg.SRe = real(specgram);
sg.SIm = imag(specgram);
sg.target = x;
identity = full(sol3(:));
func1 = @(x) corrRealMtxFastExtendedPagedHC2(x, origHalfLen, origFFTLen, newHalfLen, Ls, hop, sg, [], outterPad);
[fval, grad, y1] = func1(identity);
% plot(y1 - x);axis tight
fvtool(y1 - x)
%%
desiredFrameCount = Ls / 16;
cpLen = ceil((desiredFrameCount + (Ls / hop / 2)) * hop);
numDemoSignal = 16;
for idx = 1 : numDemoSignal
    target = randn(cpLen, 1);
    specgram = ltv_spectrogram2(target, coeff, outterPad);
    specgram(origHalfLen+1:origFFTLen, :) = conj(specgram(origHalfLen-1:-1:2, :));
    sg2.SRe(:, :, idx) = real(specgram);
    sg2.SIm(:, :, idx) = imag(specgram);
    sg2.target(:, idx) = target;
end
func2 = @(x) corrRealMtxFastExtendedPaged2(x, origHalfLen, origFFTLen, newHalfLen, Ls, hop, sg, [], outterPad);
func2_ = @(x) corrRealMtxFastExtendedPagedHC2(x, origHalfLen, origFFTLen, newHalfLen, Ls, hop, sg2, gWeights(:), outterPad);
opt.optTol = 1e-60;
opt.progTol = 1e-60;
opt.MaxIter = 1000;
opt.MaxFunEvals = opt.MaxIter * 2;
opt.Method = 'scg';
weightsOpt2 = minFunc(func2_, (identity), opt);
[fval2_, grad2_, preview2_] = func1(weightsOpt2);
weightsOpt2 = reshape(weightsOpt2, newHalfLen, origFFTLen);
% plot(y1 - x);axis tight
% fvtool(y1 - x)
%% Real time loop
xPadded = [zeros(size(Phi_s2, 1) * 0, 1); x(:); zeros(size(Phi_s2, 1), 1)];
nframes = fix(ceil(length(xPadded)/hop) - fftLen / hop / 2);
% zero padding at the end to complete the last frame
paddedLen = hop * nframes;
xPadded = [xPadded; zeros(length(xPadded) - paddedLen + hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * hop;
y = zeros(nframes * hop, 1);
accFrame2 = zeros(Ls, 1);
shtfft2 = -(2*rem(0:(origHalfLen-1), 2) - 1)';
for m = 1 : nframes % index of window position
    % Xa2 = LSsolution * fft(x(frmIdx(m) : frmIdx(m) + fftLen - 1)); % compute analysis filter outputs of a block of L samples
    Xa2 = Phi0 * x(frmIdx(m) : frmIdx(m) + fftLen - 1); % compute analysis filter outputs of a block of L samples
    % Xa2 = Xa2(1 : origHalfLen);
    % ym = Phi_s2*Xa2; % reconstruct a block of Ls samples
    Xa = fft(fftshift(xPadded(frmIdx(m) : frmIdx(m) + fftLen - 1) .* analysisWnd)); % compute analysis filter outputs of a block of L samples
    q_fft_frameinv = ltv_1st_2ndNoSq2(Xa, b, a, c1, c2);
    lastElement = q_fft_frameinv(end);
    q_fft_frameinv = q_fft_frameinv(1 : origHalfLen + 1);
    q_fft_frameinv = 0.5 * (q_fft_frameinv + 0.5 * [lastElement; q_fft_frameinv(1 : end - 1)] + 0.5 * [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
    q_fft_frameinv2 = shtfft2 .* q_fft_frameinv(1 : origHalfLen);
    q_fft_frameinv2(1) = real(q_fft_frameinv2(1));
    q_fft_frameinv2(end) = real(q_fft_frameinv2(end));
    q_fft_frameinv2(origHalfLen+1:origFFTLen) = conj(q_fft_frameinv2(origHalfLen-1:-1:2));
    plot(abs(Xa2 - q_fft_frameinv2))
    axis tight
    disp('')
    q_fft_frameinv3 = sol3 * q_fft_frameinv2;
    q_fft_frameinv3(1) = real(q_fft_frameinv3(1));
    q_fft_frameinv3(end) = real(q_fft_frameinv3(end));
    q_fft_frameinv3(halfLen2+1:size(Phi_s2, 1)) = conj(q_fft_frameinv3(halfLen2-1:-1:2));
    % Overlap add
    accFrame2 = accFrame2 + ifft(q_fft_frameinv3);
    myOut = accFrame2(1 : hop);
    accFrame2 = [accFrame2(hop + 1 : end); zeros(hop, 1)];
    y(frmIdx(m):frmIdx(m)+hop-1) = myOut;
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
coeff.origFFTLen = fftLen;
coeff.halfLen = halfLen;
coeff.hop = hop;
coeff.fs = fs;
coeff.reqSynthesisWnd = reqSynthesisWnd;
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