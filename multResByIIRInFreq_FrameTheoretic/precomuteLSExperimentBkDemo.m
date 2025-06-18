function [coeff, f] = precomuteLSIIRBkDemo(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
if nargin < 1
    fs = 48000;
    fftLen = 40;
    hop = 17;
    oct = 5;
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
prepad = fix(fftLen / 4);
pospad = fix(fftLen / 4);
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
    thetas1 = [fliplr(thetas(2 : prepad + 1)), thetas, fliplr(thetas(end - pospad + 1: end - 1))];
else
    thetas1 = [fliplr(thetas(2 : prepad + 1)), thetas, fliplr(thetas(end - pospad + 1: end - 1))];
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
theoreticalWindowShape = [zeros(size(thetas1, 2), 1), h2(:, (halfLen-1):-1:2), h2];
theoreticalWindowShape = theoreticalWindowShape(prepad + 1 : end - pospad + 1, :);
theoreticalWindowShape = [zeros(size(theoreticalWindowShape, 1), eachSide), theoreticalWindowShape, zeros(size(theoreticalWindowShape, 1), eachSide)];
%%
theoreticalWindowShape1 = theoreticalWindowShape(:, eachSide + 1 : eachSide + fftLen);
% for idx = 1 : size(theoreticalWindowShape, 1)
%     theoreticalWindowShape1(idx, :) = analysisWnd.^10;
% end
theoreticalWindowShape2 = theoreticalWindowShape1;
theoreticalWindowShape2(halfLen+1:fftLen,:) = conj(theoreticalWindowShape2(halfLen-1:-1:2,:));
centred2 = theoreticalWindowShape2;
dftMtx = dftmtx(fftLen);
windowedDFTMatrix2 = dftMtx .* centred2;
LSsolution = ifft(windowedDFTMatrix2.');
% LSsolution = real(LSsolution);
% LSsolution(abs(LSsolution) < 1e-8) = 0;
% kk = real(conj(dftMtx) .* fft(LSsolution)).';
Phi0 = LSsolution * dftMtx;
Phi0 = Phi0(1 : halfLen, :);

tol = 1e-14;
% build frame matrix Phi of Xa = Phi*x
% Phi is used to compute pseudo inverse (dual frame matrix)
ckmin = 8;
ckmin = 2;
ck = max(ckmin,floor(fftLen/halfLen)); % determines size of matrix Phi (and reconstruction error)
Nx = ck*fftLen;
Nm = floor(Nx/hop);
k = ck*floor(fftLen/hop);
ic = halfLen*floor(k*halfLen/(2*halfLen)); % index to select matrix block
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
w = size(tPhi0,2);
h = size(tPhi0,1);
N = Nm;
vshft = hop;
hshft = halfLen;
nRows = vshft*(N-1) + h;
nCols = hshft*(N-1) + w;
tPhi = sparse(nRows,nCols);
[jdx0,idx0] = meshgrid(1:w,1:h);
for l = 1:N
    idx = idx0 + (l-1)*vshft;
    jdx = jdx0 + (l-1)*hshft;
    % in case you want really overlapping to add overlapping blocks
    tPhi = tPhi + sparse(idx(:),jdx(:),tPhi0(:),nRows,nCols);
end
blk = real(Phi0'*Phi0);
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
    % in case you want really overlapping to add overlapping blocks
    C = C + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
end
N = fftLen / hop + 1;
nRows = vshft*(N-1) + h;
nCols = hshft*(N-1) + w;
pruned = sparse(nRows,nCols);
[jdx0,idx0] = meshgrid(1:w,1:h);
for l = 1:N
    idx = idx0 + (l-1)*vshft;
    jdx = jdx0 + (l-1)*hshft;
    % in case you want really overlapping to add overlapping blocks
    pruned = pruned + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
end
toc
fprintf(1,'size of full analysis frame matrix = %dx%d\n',fliplr(size(tPhi)));
fprintf(1,'computing inverse analysis frame matrix...\n');
%%
epsi = 1e-14;
% C2 = full(C);
ok = C+epsi*speye(size(C));
Phi_s = ok\tPhi;  % pseudo inverse (faster than pinv(Phi)
padUp = ic / halfLen * hop;
padDown = hop*(Nm-1) - padUp;
tPhiSel = full(tPhi(:, ic+1:ic+halfLen));
disp([padUp, fftLen])
bb = [zeros(padUp, halfLen); tPhi0; zeros(padDown, halfLen)];
disp(sum(sum(abs(tPhiSel-bb))))
pr = tPhiSel(fftLen+1:ck*fftLen, :);
padUp = ic / halfLen * hop - fftLen;
if padUp < 1
    cut1 = (ck-1)*fftLen - fftLen - padUp;
else
    cut1 = 0;
end
cut2 = (ck-1)*fftLen - fftLen - padUp;
cc = [zeros(padUp, halfLen); tPhi0(cut1 + 1 : end, :); zeros(cut2, halfLen)];
disp(sum(sum(abs(pr-cc))))
Phi_cs = ok(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen))\tPhiSel(fftLen+1:ck*fftLen, :); % pseudo inverse (faster than pinv(Phi)
R = cholesky(ok(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen)));
figure(1)
imagesc(log10(abs(R) + eps))
R2 = cholesky(pruned);
figure(2)
imagesc(log10(abs(R2 + eps)))
R3 = R2(1 : fftLen + hop * (N - 2), 1 : fftLen + hop * (N - 2));
figure(3)
imagesc(log10(abs(R3 + eps)))
% R2 = chol(ok);
% Phi_s2 = R2 \ (R2' \ tPhi);
% res = backwardSubstitution(R', forwardSubstitution(R, full(tPhi(:, 2))));
% since Phi is complex
% regularization is needed for larger tol

% perform matrix factorization of a block of Phi_s
% note: block must be chosen after transition period (L samples)
Phi_s1 = Phi_s(fftLen+1:ck*fftLen, :); % inverse matrix without transition periods
Phi_s2 = Phi_s1(:, ic+1:ic+halfLen); % matrix block which is repeated in inverse frame operator, and used in synthesis stage
Phi_s3 = Phi_cs(fftLen+1:ck*fftLen, :); % matrix block which is repeated in inverse frame operator, and used in synthesis stage
% limit Phi_s2 according to chosen reconstruction error
tol = tol*max(abs((Phi_s2(:,1))));
imin = find(abs(Phi_s2(:,1))>tol, 1,'first');
imax = find(abs(Phi_s2(:,1))>tol, 1,'last');
Phi_s2 = Phi_s2(imin:imax,:);    % limit size of reconstruction matrix
Phi_s2 = full(Phi_s2);
Ls = size(Phi_s2,1);
fprintf(1,'size of synthesis frame matrix = %dx%d\n',Ls,halfLen);
x = randn(16384, 1);
Nx = length(x);
nr = size(x,1);
x = [zeros(fftLen*0,1); x(:) ; zeros(fftLen*4,1)];
Nm = floor(Nx/hop)+1;
y = zeros(Nx+Ls,1);
nx = 1:fftLen;
ny = 1:Ls;
f = (0:1:fftLen/2)*fs/fftLen;
shtfft = -(2*rem(0:(halfLen-1), 2) - 1)';
for m = 0:Nm-1 % index of window position
    Xa2 = LSsolution * fft(x(nx+m*hop)); % compute analysis filter outputs of a block of L samples
    Xa2 = Xa2(1 : halfLen);
    Xa = fft(fftshift(x(nx+m*hop) .* analysisWnd)); % compute analysis filter outputs of a block of L samples
    Xa = Xa(1 : halfLen);
    q_fft_frameinv = shtfft .* ltv_1st_2ndNoSq2(Xa, b, a, c1, c2, prepad, pospad, [], halfLen);
    plot(real(Xa2) - real(q_fft_frameinv));
    axis tight
    % Xa2 = Xa;
    % Xa2(halfLen+1:fftLen,:) = conj(Xa2(halfLen-1:-1:2,:));
    ym = Phi_s2*Xa2; % reconstruct a block of Ls samples
    % plot(x(nx+m*hop));hold on
    % plot(real(ym));
    % hold off;axis tight
    % from subband signals Xa
    y(ny+m*hop) = y(ny+m*hop)+ym;
end

x = x(1:Nx);
if isreal(x)
    y = real(y);
end
[~,imax] = max(abs(xcorr(x,y(1:Nx))));
Nd = Nx-imax;                    % signal delay (used to compute error)
if Nd>0
    y = y(Nd+1:Nx+Nd);
else
    y = [zeros(-Nd,1) ; y(1:Nx+Nd)];
end
e = x-y;
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
function [L, r] = cholesky(A)
[~, n] = size(A);
% Full rank Cholesky factorization of A
dA = diag(A);
tol = min(dA(dA > 0)) * 1e-9;
L = zeros(size(A));
r = 0;
for k = 1 : n
    r = r + 1;
    L(k : n, r) = A(k : n, k) - L(k : n, 1 : (r - 1)) * L(k, 1 : (r - 1))';
    % Note: for r=1, the substracted vector is zero
    if (L(k, r) > tol)
        L(k, r) = sqrt(L(k, r));
        if (k < n)
            L((k + 1) : n, r) = L((k + 1) : n, r) / L(k, r);
        end
    else
        r = r - 1;
    end
end
end
function x = forwardSubstitution(L, b)
%This function solve a lower triangular system using forward substitution
%method. The standard call is: "x=forwardSubstitution(U,b) in wich L and b represent
%respectively the lower triangular matrix and  the known term.
%******************
%Riccardo Dessì
%e-mail:ri.dessi1@studenti.unica.it
%******************
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
%This function solve an up triangular system using backwardSubstitution substitution
%method. The standard call is: "x=backwardSubstitution(U,b) in wich U and b represent
%respectively the up triangular matrix and  the known term.
%******************
%Riccardo Dessì
%e-mail:ri.dessi1@studenti.unica.it
%******************
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