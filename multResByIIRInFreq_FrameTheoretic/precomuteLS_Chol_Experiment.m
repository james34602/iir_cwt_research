function [coeff, f] = precomuteLS_Chol(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 128;
    hop = 27;
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
for idx = 1 : size(theoreticalWindowShape1, 1)
    Xa2 = circshift(Xa, idx - 1);
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
 Phi0(:)=randn(size(Phi0));
Phi0(1, :) = real(Phi0(1, :));
Phi0(end, :) = real(Phi0(end, :));
Phi0(halfLen+1:fftLen, :) = conj(Phi0(halfLen-1:-1:2, :));
origHalfLen = halfLen;
origFFTLen = fftLen;
[halfLen, fftLen] = size(Phi0);
% build frame matrix Phi of Xa = Phi*x
% Phi is used to compute pseudo inverse (dual frame matrix)
ckmin = 4;
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
tPhi0 = Phi0';
blk = real(tPhi0*Phi0);
w = size(blk,2);
h = size(blk,1);
N = Nm + 1;
vshft = hop;
hshft = hop;
nRows = vshft*(N-1) + h;
nCols = hshft*(N-1) + w;
C = sparse(nRows,nCols);
[jdx0,idx0] = meshgrid(1:w,1:h);
for l = 1:N
    idx = idx0 + (l-1)*vshft;
    jdx = jdx0 + (l-1)*hshft;
    if l > 1
        C = C + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    else
        C = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    end
end
% C1 = zeros(hop*(N-1) + fftLen - hop, hop*(N-1) + fftLen - hop);
% for l = 1:N-1
%     if l > 1
%         C1((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) = C1((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) + blk;
%     else
%         C1((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) = blk;
%     end
% end
% C3 = C1(hop*floor((fftLen-offset)/hop)-hop + 1: hop*floor((fftLen-offset)/hop)-hop + fftLen, hop*floor((fftLen-offset)/hop)-hop + 1: hop*floor((fftLen-offset)/hop)-hop + fftLen);
C1_ = zeros(fftLen, fftLen);
validRegion = fftLen - (hop*floor((fftLen-offset)/hop)-hop);
for l = 1:N-1
    if l > 1
        if validRegion <= fftLen
            dest = 1 : validRegion;
            selBlk = blk((fftLen - validRegion) + 1 : fftLen, (fftLen - validRegion) + 1 : fftLen);
        else
            dest = (validRegion - fftLen) + 1 : fftLen;
            selBlk = blk(1 : (fftLen - (validRegion - fftLen)), 1 : (fftLen - (validRegion - fftLen)));
        end
        C1_(dest, dest) = C1_(dest, dest) + selBlk;
    else
        C1_(1 : validRegion, 1 : validRegion) = blk((fftLen - validRegion) + 1 : fftLen, (fftLen - validRegion) + 1 : fftLen);
    end
    validRegion = validRegion + hop;
end
% sum(abs(C1_(:)-C3(:)))
C3 = C1_;
epsi = 1e-14;
C=full(C);
Cminimal = C(fftLen-offset+1:ck*fftLen-offset, (fftLen-offset+1:ck*fftLen-offset)) + epsi*speye(ck*fftLen-fftLen);
C2 = C(hop*floor((fftLen-offset)/hop) + 1: hop*ceil((ck*fftLen-offset)/hop), hop*floor((fftLen-offset)/hop) + 1: hop*ceil((ck*fftLen-offset)/hop));
solCholSize = hop*ceil((ck*fftLen-offset)/hop) - hop*floor((fftLen-offset)/hop);
% return;
C2 = C2 + epsi*speye(size(C2));
C3 = C3 + epsi*speye(size(C3));
pd = fftLen-offset - (hop*floor((fftLen-offset)/hop));
% imagesc(log(abs(C2)))
padUp = ceil(k*halfLen/(2*halfLen)) * hop + offset;
padDown = hop*(Nm-1) - padUp;
tPhiSel = [zeros(padUp, halfLen); tPhi0; zeros(padDown, halfLen)];
%%
R1 = chol(C2, 'upper');
R2_ = chol_sec_nr_blockShortCut3(C3, solCholSize, fftLen, hop);
% R2 = chol_sec_nr_blockShortCut2(full(C2), fftLen, hop);
% R2 = chol_sec_nr_blockShortCut(C6, hop, fftLen, hop);
R3 = chol_sec_nr_block(C2, hop, fftLen, hop);
R4 = croutChol(C2)';
%% Cholesky
cuttedPos = tPhiSel(fftLen+1:ck*fftLen, :);
Phi_s2 = Cminimal \ cuttedPos; % pseudo inverse
trim = ((hop*ceil((ck*fftLen-offset)/hop)) - (hop*floor((fftLen-offset)/hop))) - (((ck*fftLen)-fftLen) + pd);
cuttedPos2 = [zeros(pd, halfLen); cuttedPos; zeros(trim, halfLen)];
Phi_s4 = R3 \ (R3' \ cuttedPos2); % pseudo inverse
Phi_s3 = Phi_s4(pd + 1 : pd + ((ck*fftLen)-fftLen), :);
Phi_s2 = Phi_s3;

Ls = size(Phi_s2,1);
newHalfLen = getFFTHalfLen(size(Phi_s2, 1));
% dmt = dftmtx(size(Phi_s2, 1))';
% sol2 = real(dmt\Phi_s2 * size(Phi_s2, 1));
sol2 = real(fft(Phi_s2, [], 1));
halfLen2 = size(Phi_s2, 1) / 2 + 1;
sol3 = sol2(1 : halfLen2, :);

sigLen = 32768;
x = randn(sigLen, 1);
%%
coeff.origHalfLen = origHalfLen;
coeff.origFFTLen = origFFTLen;
coeff.Ls = Ls;
coeff.hop = hop;
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
end
function U = chol_sec_nr_blockShortCut3(submatrix, solCholSize, fftLen, hop)
U = zeros(solCholSize, solCholSize);
inds = 1:hop;
% Loop through each chunk of blocksize:
submtxBk = submatrix;
remaining = solCholSize;
opts.LT = true;
for n = 1 : solCholSize / hop
    remaining = remaining - hop;
    myOffset = (n-1)*hop;
    Utl = chol(submatrix(inds, inds));
    if remaining < fftLen - hop
        validLen = remaining;
    else
        validLen = fftLen - hop;
    end
    idxorig = hop + 1 : hop + validLen;
    idxSbr = myOffset+hop + 1 : myOffset+hop + validLen;
    selectedA = submatrix(inds, idxorig);
    Utr = linsolve(Utl', selectedA, opts); %Utl' \ selectedA;

    U(inds + myOffset,inds+myOffset) = Utl;
    U(myOffset + inds, idxSbr) = Utr;

    gram = Utr' * Utr;
    % imagesc(log(abs(U)))
    % clim([-3, 5.2])
    idxRelative = 1 : validLen;
    idxRelative2 = 1 : (validLen + hop);
    % imagesc(log(abs(submatrix)))
    nsubmtx = submatrix(hop + idxRelative, hop + idxRelative) - gram;
    submatrix = submtxBk(idxRelative2, idxRelative2);
    submatrix(idxRelative, idxRelative) = nsubmtx;
end
end
function U = chol_sec_nr_blockShortCut2(A, fftLen, hop)
if(isempty(A))
    U = [];
    return;
end
N = length(A);
U = zeros(size(A));
inds = 1:hop;
% Loop through each chunk of blocksize:
submatrix = A(1 : fftLen, 1 : fftLen);
submtxBk = submatrix;
remaining = N;
for n = 1:N/hop
    remaining = remaining - hop;
    myOffset = (n-1)*hop;
    Utl = chol(submatrix(inds, inds));
    if remaining < fftLen - hop
        validLen = remaining;
    else
        validLen = fftLen - hop;
    end
    idxorig = hop + 1 : hop + validLen;
    idxSbr = myOffset+hop + 1 : myOffset+hop + validLen;
    selectedA = submatrix(inds, idxorig);
    Utr = Utl' \ selectedA;

    U(inds + myOffset,inds+myOffset) = Utl;
    U(myOffset + inds, idxSbr) = Utr;

    gram = Utr' * Utr;
    % imagesc(log(abs(U)))
    % clim([-3, 5.2])
    idxRelative = 1 : validLen;
    idxRelative2 = 1 : (validLen + hop);
    % imagesc(log(abs(submatrix)))
    nsubmtx = submatrix(hop + idxRelative, hop + idxRelative) - gram;
    submatrix = submtxBk(idxRelative2, idxRelative2);
    submatrix(idxRelative, idxRelative) = nsubmtx;
end
end
function U = chol_sec_nr_blockShortCut(A, blockSize, fftLen, hop)
Aorig=A;
if(isempty(A))
    U = [];
    return;
end
N = length(A);
U = zeros(size(A));
inds = 1:blockSize;
% Loop through each chunk of blocksize:
submatrix = A(blockSize + 1 : blockSize + fftLen, blockSize + 1 : blockSize + fftLen);
submtxBk = submatrix;
remaining = N;
for n = 1:N/blockSize
    remaining = remaining - blockSize;
    myOffset = (n-1)*blockSize;

    % imagesc(log(abs(inp)))
    curIdx = myOffset + inds;
    nextIdx = n*blockSize + inds;
    inp = A(myOffset + inds, myOffset + inds);
    inp2 = submatrix(inds, inds);
    disp(sum(abs(inp(:) - inp2(:))))
    Utl = chol(A(myOffset + inds, myOffset + inds));
    % Utl2 = chol(submatrix(inds, inds));
    if remaining < fftLen - hop
        validLen = remaining;
    else
        validLen = fftLen - hop;
    end
    if (remaining - blockSize) < fftLen - hop
        validLen2 = remaining - blockSize;
    else
        validLen2 = fftLen - hop;
    end
    idxorig = blockSize + 1 : blockSize + validLen;
    idxSbr = myOffset+blockSize + 1 : myOffset+blockSize + validLen;
    idxnext = (n)*blockSize+blockSize + 1 : (n)*blockSize+blockSize + validLen2;
    idxRelative = 1 : validLen;
    selectedA = A(myOffset + inds, idxSbr);
    selectedA_ = submatrix(inds, idxorig);
    disp(sum(abs(selectedA_(:) - selectedA(:))))
    Utr = Utl' \ selectedA;

    U(inds + myOffset,inds+myOffset) = Utl;
    U(myOffset + inds, idxSbr) = Utr;

    gram = Utr' * Utr;
    % imagesc(log(abs(A(idxSbr, idxSbr))))
    % clim([-3, 5.2])
    Abk = A(idxSbr, idxSbr);
    disp(sum(sum(abs(A(idxSbr, idxSbr) - submatrix(blockSize + idxRelative, blockSize + idxRelative)))))
    if n == 8
        disp('')
    end
    A(idxSbr, idxSbr) = A(idxSbr, idxSbr) - gram;
    % Anext_ = Aorig(idxnext, idxnext);
    % Anext = A(idxnext, idxnext);
    % if remaining < fftLen
    %     validLen3 = fftLen - (fftLen - remaining);
    % else
    %     validLen3 = fftLen;
    % end
    % submtxBk = Aorig(n * blockSize + 1 : n * blockSize + validLen3, n * blockSize + 1 : n * blockSize + validLen3);
    % imagesc(log(abs(submtxBk)))
    % imagesc(Anext_ - Anext)
    gt = A(idxSbr, idxSbr);
    if n == 7
        disp('')
    end
    nsubmtx = submatrix(blockSize + idxRelative, blockSize + idxRelative) - gram;
    submatrix = submtxBk;
    submatrix(idxRelative, idxRelative) = nsubmtx;
    % imagesc(log(abs(A(myOffset+1+blockSize:N, myOffset+1+blockSize:N))))
    if n == 7
        disp('')
    end
end
end
function U = chol_sec_nr_block(A, blockSize, fftLen, hop)
Abk = A;
if(isempty(A))
    U = [];
    return;
end
N = length(A);
U = zeros(size(A));

inds = 1:blockSize;

% Loop through each chunk of blocksize:
lst1 = [];
histInp = zeros(blockSize, blockSize, N/blockSize);
remaining = N;
for n = 1:N/blockSize
    remaining = remaining - blockSize;
    myOffset = (n-1)*blockSize;
    inp = A(myOffset + inds, myOffset + inds);
    histInp(:, :, n) = inp;
    lst1 = [lst1; myOffset + inds];

    % aa = croutChol(A(myOffset + inds, myOffset + inds));
    % imagesc(log(abs(inp)))
    Utl = chol(A(myOffset + inds, myOffset + inds));
    if remaining < fftLen - hop
        validLen = remaining;
    else
        validLen = fftLen - hop;
    end
    idxSbr = myOffset+blockSize + 1 : myOffset+blockSize + validLen;
    selectedA = A(myOffset + inds, idxSbr);
    Utr = Utl' \ selectedA;

    U(inds + myOffset,inds+myOffset) = Utl;
    U(myOffset + inds, idxSbr) = Utr;

    gram = Utr' * Utr;
    % imagesc(log(abs(A(idxSbr, idxSbr))))
    A(idxSbr, idxSbr) = A(idxSbr, idxSbr) - gram;
    % imagesc(log(abs(A(myOffset+1+blockSize:N, myOffset+1+blockSize:N))))
end
end
function U = chol_sec_nr_blockOrig(A, blockSize)
if(isempty(A))
    U = [];
    return;
end
N = length(A);
U = zeros(size(A));

myOffset = -blockSize;
inds = 1:blockSize;

% Loop through each chunk of blocksize:
for n = 1:N/blockSize
    myOffset = (n-1)*blockSize;

    Utl = chol(A(myOffset + inds, myOffset + inds));
    Utr = Utl'\A(myOffset + inds, (n*blockSize+1):N);

    U(inds + myOffset,inds+myOffset) = Utl;
    U(myOffset + inds, (n*blockSize+1):N) = Utr;

    Sbr = A(myOffset+1+blockSize:N, myOffset+1+blockSize:N) - Utr'*Utr;
    % imagesc(log(abs(U)))
    A(myOffset+1+blockSize:N, myOffset+1+blockSize:N)  = Sbr;
end
% finish off the rest of the computation < blockSize int he bottom right:
if(mod(N, blockSize) > 0)
    U(myOffset+1+blockSize:end, myOffset+1+blockSize:end) = chol(A(myOffset+1+blockSize:end, myOffset+1+blockSize:end));
end
end
function L = croutChol(A)
[~, n] = size(A);
L = zeros(size(A));
for j = 1 : n
    diagonal_value = sqrt(A(j, j) - L(j, j));
    L(j, j) = diagonal_value;

    diagonal_value = 1.0 / diagonal_value;

    % #pragma omp parallel for private(i,k,sum) shared (A,L,j,diagonal_value) schedule(static) if (j < n - cutoff)
    for i = j + 1 : n
        sum = L(i, 1 : j) * L(j, 1 : j)';
        value = (diagonal_value * (A(i, j) - sum));
        L(i, j) = value;
        L(i, i) = L(i, i) + value * value;
    end
end
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
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