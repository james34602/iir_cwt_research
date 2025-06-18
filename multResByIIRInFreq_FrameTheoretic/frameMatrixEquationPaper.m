function frameMatrixEquationPaper(N, L, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    N = 64;
    L = 8;
    oct = 20;
    order = 2;
    reqSynthesisWnd = 1;
    sparCutOff = 5e-8;
    zp = 1;
end
outterWnd = 1;
paddedFFTLen = N * zp;
eachSide = (paddedFFTLen - N) / 2;
M = getFFTHalfLen(N);
paddedHalfLen = getFFTHalfLen(paddedFFTLen);
ovp = N / L;
% number of points of pre and post padding used to set initial conditions
prepad = fix(N / 2);
pospad = fix(N / 2);
% frequency bins grid (linear in this case) - pre and pos padding is added
% poles of the IIR LTV Q FFT transform for the parameters above
f = (0:1:N/2)*fs/N;
% number of points of pre and post padding used to set initial conditions
thetas = 0:(N/2);
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
analysisWnd = hann(N, 'periodic') .^ (1 / reqSynthesisWnd);
synWnd1 = hann(N, 'periodic') .^ reqSynthesisWnd;
chopedWnd1 = analysisWnd(M : end);
chopedWnd2 = synWnd1(M : end);
halfWndLen = M - 1;
digw2 = linspace(0, pi, M);
digw = digw2(1 : halfWndLen);
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd1.';
h2 = h2 .* chopedWnd2.';
theoreticalWindowShape = [zeros(size(thetas, 2), 1), h2(:, (M-1):-1:2), h2];
tol = min(1.5, max(1, (N / L) / 5));
hopsizeTol = min(N, ceil(L * 2 * tol));
if mod(hopsizeTol, 2) == 1
    hopsizeTol = hopsizeTol - 1;
end
smallestPossibleWnd = [zeros((N - hopsizeTol) / 2, 1); hann(hopsizeTol, 'periodic'); zeros((N - hopsizeTol) / 2, 1)];
wndDif = theoreticalWindowShape' - smallestPossibleWnd;
wndDifPwr = sum(abs(wndDif), 1);
[~, firstUndersampling1] = min(wndDifPwr);
firstUndersampling2 = find(any(wndDif < 0), 1, 'first');
firstUndersampling = ceil((firstUndersampling1 + firstUndersampling2) / 2);
firstUndersampling = max(firstUndersampling, fix(N / L * oct / 2));
if ~isempty(firstUndersampling)
    thetas = 0:(N/2);
    thetas(1) = eps;
    thetas(firstUndersampling : end) = thetas(firstUndersampling);
    thetas1 = thetas;
    thetas1(M+1:N) = conj(thetas1(M-1:-1:2));
else
    thetas1 = thetas;
    thetas1(M+1:N) = conj(thetas1(M-1:-1:2));
end
% Eliminate oscillation around corner
time = 0.026 * N; % More elements in array less smoothing is needed
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
sigmas = (thetas1 ./ N) ./ oct / pi * N;
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
%% Obtain DFT filterbank frequency response
if order == 2
    h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
else
    h = (cplxFreq .* a(:)) ./ (cplxFreq - (1 - a(:)));
end
h2 = (h .* conj(h)) .* chopedWnd1.';
h2 = h2 .* chopedWnd2.';
theoreticalWindowShape = [zeros(size(thetas1, 2), 1), h2(:, (M-1):-1:2), h2];
theoreticalWindowShape = [zeros(size(theoreticalWindowShape, 1), eachSide), theoreticalWindowShape, zeros(size(theoreticalWindowShape, 1), eachSide)];
%%
theoreticalWindowShape1 = theoreticalWindowShape(:, eachSide + 1 : eachSide + N);
actualWindowShape3 = zeros(size(theoreticalWindowShape1));
dftMtx = dftmtx(N);
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
shtfft = -(2*rem(0:(M-1), 2) - 1)';
for idx = 1 : size(theoreticalWindowShape1, 1)
    Xa2 = circshift(Xa, idx - 1);
    % Xa2_ = Xa3(:, idx);
    % Xa2 and Xa2_ is equal symbolically, but not numerically exact
    % disp(sum(abs(Xa2 - Xa2_)))
    yy = ltv_1st_2ndNoSq2(Xa2, b, a, c1, c2, prepad, pospad, [], M);
    b2 = fftshift(conj(dftMtx(:, idx)) .* fft(yy));
    % plot(real(b2));
    % hold on
    % plot(imag(b2))
    % hold off
    % axis tight
    actualWindowShape3(idx, :) = b2;
end
actualWindowShape3 = (actualWindowShape3 .* synWnd1.') / N;
windowedDFTMatrix2 = dftMtx .* actualWindowShape3;
LSsolution = ifft(windowedDFTMatrix2.');
% LSsolution = real(LSsolution);
dmt = dftmtx(320);
% kk = real(conj(dftMtx) .* fft(LSsolution)).';
eachSide = (320 - N) / 2;
Phi0 = LSsolution * dftMtx;
% Phi1 = [zeros(fftLen, eachSide), Phi0, zeros(fftLen, eachSide)]';
% gg = dmt \ Phi1;
% Phi0 = Phi1;
Phi0 = Phi0(1 : M, :);
Phi0(1, :) = real(Phi0(1, :));
Phi0(end, :) = real(Phi0(end, :));
Phi0(M+1:N, :) = conj(Phi0(M-1:-1:2, :));
LSsolution = real(Phi0 / dftMtx);
origHalfLen = M;
origFFTLen = N;
% halfLen = fftLen;
[M, N] = size(Phi0);

tol = 1e-16;
% build frame matrix Phi of Xa = Phi*x
% Phi is used to compute pseudo inverse (dual frame matrix)
ckmin = 6;
% ckmin = 2;
ck = max(ckmin,floor(N/M)); % determines size of matrix Phi (and reconstruction error)
Nx = ck*N;
Nl = floor(Nx/L);
k = ck*floor(N/L);
ic = M*floor(k*M/(2*M)); % index to select matrix block
if L*(Nl-1) + N < ck*N
    disp([L*(Nl-1) + N, ck*N])
    error('Synthesis frame matrix does''t exist');
end
% tic
Phi_ = zeros(M*Nl,(Nl-1)*L+M);
for m = 0:Nl-1
    Phi_(m*M+1:m*M+M,m*L+1:m*L+M) = Phi0;
end
S = zeros(M*Nl,(Nl-1)*L+N);
D = kron(sparse(eye(Nl)), Phi0);
I = eye(M, N);
for m = 0:Nl-1
    i = m*M + (1 : M);
    j = m*L + (1 : N);
    S(i, j) = I;
end
S2 = shifted_identity_matrix(N, M, L, Nl);
Phi = D * S;
end
function S = shifted_identity_matrix(M, N, L, Nl)
% Initialize the final matrix S with zeros
total_rows = Nl * M;  % Number of rows
total_cols = (Nl - 1) * L + N;  % Number of columns
S = zeros(total_rows, total_cols);

% Loop over each block placement (m = 0, 1, ..., Nl-1)
for m = 0:(Nl-1)
    % Create E_m(i,j) as identity matrix I_{M x N}
    E_m = eye(M, N);  % Identity block

    % Define the row and column indices for placing the identity matrix
    row_start = m * M + 1;         % Row index starts at m * M + 1
    row_end = row_start + M - 1;   % Row index ends at m * M + M

    col_start = m * L + 1;         % Column index starts at m * L + 1
    col_end = col_start + N - 1;   % Column index ends at m * L + N

    % Insert E_m(i,j) into S at the specified position
    S(row_start:row_end, col_start:col_end) = E_m;
end
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