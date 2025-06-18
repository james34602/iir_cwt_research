function create4ndLTVWaveletBasis(fftLen, hop, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
addpath('../')
addpath('../minFunc_2012/autoDif')
addpath('../minFunc_2012/minFunc')
addpath('../minFunc_2012/minFunc/compiled')
if nargin < 1
    fs = 48000;
    fftLen = 1024;
    hop = 16;
    oct = 2;
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
    [~, ~, ~, Bmfwd, Bmbwd, Am] = InitDeriche(sigmas(1), length(sigmas));
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
    % yy = ltv_1st_2ndNoSq2(Xa2, b, a, c1, c2);
    yy = ltv_1st_2ndNoSq1(Xa2, Bmfwd, Bmbwd, Am, divareaSum);
    b2 = fftshift(conj(dftMtx(:, idx)) .* fft(yy));
    actualWindowShape3(idx, :) = b2;
end
actualWindowShape3 = (actualWindowShape3 .* synWnd1.') / fftLen;
%% Plot result after overlap-add
clear cpxRes dftFilterbank2 theoreticalWindowShape2 theoreticalWindowShape
windowedDFTMatrix1 = dftMtx .* actualWindowShape3;
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
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