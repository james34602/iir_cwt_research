function [coeff, f] = precomuteLS2(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
if nargin < 1
    fs = 48000;
    fftLen = 128;
    hop = 16;
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
firstUndersampling = max(firstUndersampling, fix(fftLen / hop * oct / 2));
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
    LSsolution = Phi0 / dftMtx;
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
    ic = halfLen*ceil(k*halfLen/(2*halfLen)); % index to select matrix block
    if hop*(Nm-1) + fftLen < ck*fftLen
        disp([hop*(Nm-1) + fftLen, ck*fftLen])
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
    for l = 2:N
        idx = idx0 + (l-1)*vshft;
        jdx = jdx0 + (l-1)*hshft;
        % in case you want really overlapping to add overlapping blocks
        if l > 2
            C = C + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        else
            C = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        end
    end
    % N = ceil((fftLen+hop+fftLen)/hop);
    % nRows = hop*(N-1) + fftLen;
    % nCols = hop*(N-1) + fftLen;
    % pruned = sparse(nRows,nCols);
    % [jdx0,idx0] = meshgrid(1:fftLen,1:fftLen);
    % for l = 2:N
    %     idx = idx0 + (l-1)*hop;
    %     jdx = jdx0 + (l-1)*hop;
    %     % in case you want really overlapping to add overlapping blocks
    %     if l > 3
    %         pruned = pruned + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    %     else
    %         pruned = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    %     end
    % end
    epsi = 1e-14;
    % pruned = full(pruned(fftLen+hop +1:fftLen+hop+fftLen, fftLen+hop +1:fftLen+hop+fftLen)) + epsi*speye((fftLen+hop+fftLen)-(fftLen+hop));
    % effectiveLen = zeros(fftLen, 1);
    % prunedShifted = zeros(size(pruned));
    % for i = 1 : fftLen
    %     effectiveLen(i) = fftLen - i + 1;
    %     prunedShifted(:, i) = circshift(pruned(:, i), -(i - 1));
    %     prunedShifted(fftLen - i + 2 : end, i) = 0;
    % end
    C2 = C(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen)) + epsi*speye(ck*fftLen-fftLen);
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
    padUp = ic / halfLen * hop;
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
    Phi_s2 = C2 \ tPhiSel(fftLen+1:ck*fftLen, :); % pseudo inverse
    cuttedPos = tPhiSel(fftLen+1:ck*fftLen, :);
    % imagesc(abs(cuttedPos))
    [~, mostMaxIdx] = max(cuttedPos(:, 1), [], 1);
    newHalfLen = getFFTHalfLen(size(cuttedPos, 1));
    offset = newHalfLen - mostMaxIdx;
    if offset
        tPhiSel = circshift(tPhiSel, offset, 1);
        cuttedPos = tPhiSel(fftLen+1:ck*fftLen, :);
        % imagesc(abs(cuttedPos))
        offset=-offset;
        C2 = C(fftLen+offset+1:ck*fftLen+offset, (fftLen+offset+1:ck*fftLen+offset)) + epsi*speye(ck*fftLen-fftLen);
        Phi_s2 = C2 \ cuttedPos; % pseudo inverse
        disp('')
    end
    clear C C2 tPhiSel
else
    load('matlab2.mat');
end
Ls = size(Phi_s2,1);
% fresp = fft(Phi_s2, [], 1);
% aa = fft(conj(dftMtx), [], 1);
% nDftMtx = dftmtx(Ls);
% nv = conj(nDftMtx) .* fresp;
fprintf(1,'size of synthesis frame matrix = %dx%d\n',Ls,halfLen);
sigLen = 32768;
x = randn(sigLen, 1);
x = [zeros(size(Phi_s2, 1), 1); x(:); zeros(size(Phi_s2, 1), 1)];
nframes = fix(ceil(length(x)/hop) - fftLen / hop / 2);
% zero padding at the end to complete the last frame
paddedLen = hop * nframes;
x = [x; zeros(length(x) - paddedLen + hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * hop;
y = zeros(nframes * hop, 1);
accFrame2 = zeros(Ls, 1);
% shtfft = -(2*rem(0:(halfLen-1), 2) - 1)';
for m = 1 : nframes % index of window position
    Xa2 = LSsolution * fft(x(frmIdx(m) : frmIdx(m) + fftLen - 1)); % compute analysis filter outputs of a block of L samples
    Xa2 = Xa2(1 : halfLen);
    % Xa = fft(fftshift(x(nx+m*hop) .* analysisWnd)); % compute analysis filter outputs of a block of L samples
    % Xa = Xa(1 : halfLen);
    % q_fft_frameinv = shtfft .* ltv_1st_2ndNoSq2(Xa, b, a, c1, c2, prepad, pospad, [], halfLen);
    % plot(real(Xa2) - real(q_fft_frameinv));
    % axis tight
    % Xa2 = Xa;
    % Xa2(halfLen+1:fftLen,:) = conj(Xa2(halfLen-1:-1:2,:));
    ym = Phi_s2*Xa2; % reconstruct a block of Ls samples
    % plot(x(nx+m*hop));hold on
    % plot(real(ym));
    % hold off;axis tight
    % from subband signals Xa
    % Overlap add
    accFrame2 = accFrame2 + ym;
    myOut = accFrame2(1 : hop);
    accFrame2 = [accFrame2(hop + 1 : end); zeros(hop, 1)];
    y(frmIdx(m):frmIdx(m)+hop-1) = myOut;
end

x = x(size(Phi_s2, 1)+1:end);
y = y(size(Phi_s2, 1)+1:end);
x = x(1 : sigLen);
if isreal(x)
    y = real(y);
end
y = circshift(y, -64);
plot(x);hold on;plot(y);hold off;axis tight
y = y(1 : sigLen);
e = x-y;
%% Save coefficient
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
function [L, k] = choleskyComputationMap(A)
[~, n] = size(A);
% Full rank Cholesky factorization of A
dA = diag(A);
tol = min(dA(dA > 0)) * 1e-9;
L = zeros(size(A));
inMap = zeros(size(A), 'logical');
outMap = zeros(size(A), 'uint16');
for k = 1 : n
    inMap(k : n, k) = 1;
    L(k : n, k) = A(k : n, k) - L(k : n, 1 : (k - 1)) * L(k, 1 : (k - 1))';
imagesc(outMap)
    outMap(k : n, 1 : (k - 1)) = outMap(k : n, 1 : (k - 1)) + 1;
    outMap(k, 1 : (k - 1)) = outMap(k, 1 : (k - 1)) + 1;
    imagesc(outMap)
    outMap(k : n, k) = outMap(k : n, k) + 1;
    outMap(k, k) = outMap(k, k) + 1;
    L(k, k) = sqrt(L(k, k));
    if (k < n)
        L((k + 1) : n, k) = L((k + 1) : n, k) / L(k, k);
        outMap((k + 1) : n, k) = outMap((k + 1) : n, k) + 1;
    end
end
end
function [L, k] = cholesky(A, effectiveLen, fftLen, hop)
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
function x = forwardSubstitution(L, b)
S = size(L);
n = S(1);
if S(1) ~= S(2)
    error('matrix mast be square')
end
x = zeros(1, n);
if L(1, 1) ~= 0
    x(1, 1) = b(1) ./ L(1, 1);
else
    x(1, 1) = 0;
end
%forward substitution
for k=2:n
    if L(k, k)
        x1 = 1 / L(k, k) .* (b(k) - sum(L(k, k - 1 : -1 : 1) .* x(k - 1 : -1 : 1)));
    else
        x1 = 0;
    end
    x(1, k) = x1;
end
x = x';
end
function x = backwardSubstitution(U, b)
S = size(U);
n = S(1);
if S(1) ~= S(2)
    error('matrix mast be square')
end
x = zeros(1, n);
if U(n, n) ~= 0
    x(1, n) = b(end) ./ U(n, n);
else
    x(1, n) = 0;
end
%bacward substitution
for k = n - 1 : -1 : 1
    if U(k, k) ~= 0
        x1 = 1 / U(k, k) .* (b(k) - sum(U(k, k + 1 : end) .* x(k + 1 : end)));
    else
        x1 = 0;
    end
    x(k) = x1;
end
x = x';
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end