idx = 12;
longChirp = chirp2(t, D, f1, f2)';
tg = circshift(longChirp, fix((idx - 1) * (hop / numDemoSignal)));
[tg, fs] = loadSignal(7, fftLen);
[S1, t_q] = ltv_spectrogram2(tg, coeff);
S1(halfLen+1:fftLen,:) = conj(S1(halfLen-1:-1:2,:));
STime1 = ifft(S1) .* fftshift(synWnd1);
% ak = fftshift([zeros(eachSide, size(spec1, 2)); ifftshift(STime1); zeros(eachSide, size(spec1, 2))]);
STime1 = [STime1(1 : halfLen-1, :); zeros(eachSide * 2, size(S1, 2)); STime1(halfLen : end, :)];
S1 = fft(STime1);
getbackCorrectedToSpectrum1 = fft(STime1 .* correctionWndHF);
S1 = S1(1 : paddedHalfLen, :) .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum1(1 : paddedHalfLen, :) .* wndCorrectionWeightingHF;
% S1 = S1(1 : paddedHalfLen, :);
% Optimize real matrix with extended spectral(Extended spectral simulate circular convolution)
ext = [conj(S1(prepad + 1 : -1 : 2, :)); S1; conj(S1(paddedHalfLen - 1 : -1 : (paddedHalfLen - pospad + 1), :))];
clear sg
sg.S = ext;
sg.target = target;
func3_ = @(x) corrCplxMtxFastExtendedPagedHC2(x, paddedHalfLen, paddedFFTLen, prepad, pospad, hop, sg, [], synWnd2, fftLen, eachSide);
preview2_ = func3_(weightsOpt3);
function [S, t] = ltv_spectrogram2(x, coeff)
x = x(:);
x = [zeros(coeff.fftLen - coeff.hop, 1); x; zeros(coeff.fftLen / 2, 1)];
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.fftLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
x = [x; zeros(length(x) - paddedLen + coeff.hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
S = zeros(nframes, coeff.fftLen/2+1, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.fftLen - 1);
    X = fft(fftshift(frame .* coeff.analysisWnd));
    S(i, :) = ltv_1st_2ndNoSq2(X(1 : coeff.halfLen), coeff.b, coeff.a, coeff.c1, coeff.c2, coeff.prepad, coeff.pospad, coeff.corrF, coeff.halfLen);
end
t = (frmIdx - 1 + coeff.fftLen / 2) / coeff.fs;
S = permute(S, [2, 1, 3]);
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
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function w = tukeywin_periodic(N, r)
if r <= 0
    w = ones(N,1);
elseif r >= 1
    w = hann(N, 'periodic');
else
    L = N + 1;
    t = linspace(0,1,L)';
    % Defines period of the taper as 1/2 period of a sine wave.
    per = r/2; 
    tl = floor(per*(L-1))+1;
    th = L-tl+1;
    % Window is defined in three sections: taper, constant, taper
    w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th : end - 1) - 1 + per)))/2)];
end
end
function y = corrCplxMtxFastExtendedPagedHC2(weights, paddedHalfLen, paddedFFTLen, prepad, pospad, hop, sigs, gWeights, synWnd, fftLen, eachSide) % Same as function corrRealMtxFastExtendedPaged
x = reshape(weights, paddedHalfLen, (paddedHalfLen + prepad + pospad - 1) * 2);
x = complex(x(:, 1 : (paddedHalfLen + prepad + pospad - 1)), x(:, (paddedHalfLen + prepad + pospad - 1) + 1 : (paddedHalfLen + prepad + pospad - 1) * 2));
q_fft_frameinv2 = pagemtimes(x, sigs.S);
q_fft_frameinv2(1, :, :) = real(q_fft_frameinv2(1, :, :));
q_fft_frameinv2(end, :, :) = real(q_fft_frameinv2(end, :, :));
if mod(paddedFFTLen,2) == 0
    fg = 0;
    q_fft_frameinv2(paddedHalfLen+1:paddedFFTLen, :, :) = conj(q_fft_frameinv2(paddedHalfLen-1:-1:2, :, :));
else
    fg = 1;
    q_fft_frameinv2(paddedHalfLen+1:paddedFFTLen, :, :) = conj(q_fft_frameinv2(paddedHalfLen:-1:2, :, :));
end
if ~isempty(synWnd)
    mtx = circshiftCustom(ifft(q_fft_frameinv2), paddedFFTLen / 2) .* synWnd;
else
    mtx = circshiftCustom(ifft(q_fft_frameinv2), paddedFFTLen / 2);
end
y3 = overlapAdd3D(mtx, hop);
y = y3(fftLen - hop + eachSide + 1 : end, :);
end
function b = circshiftCustom(a,p)
m = size(a,1);
b = a(mod((0:m-1)-double(rem(p,m)), m)+1, :, :);
end