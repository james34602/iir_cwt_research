function create2ndLTVWaveletBasis(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 512;
    hop = 1;
    oct = 1;
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
windowedDFTMatrix1 = dftMtx .* actualWindowShape3; % Basis
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end