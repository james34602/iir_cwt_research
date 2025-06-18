idx = 12;
longChirp = chirp2(t, D, f1, f2)';
tg = circshift(longChirp, fix((idx - 1) * (hop / numDemoSignal)));
[tg, fs] = loadSignal(7, fftLen);
[S1, t_q] = ltv_spectrogram3(tg, coeff, LSsolution);
S1(halfLen+1:fftLen,:) = conj(S1(halfLen-1:-1:2,:));
STime1 = ifft(S1);
% ak = fftshift([zeros(eachSide, size(spec1, 2)); ifftshift(STime1); zeros(eachSide, size(spec1, 2))]);
STime1 = [STime1(1 : halfLen-1, :); zeros(eachSide * 2, size(S1, 2)); STime1(halfLen : end, :)];
S1 = fft(STime1);
% S1 = S1(1 : paddedHalfLen, :);
% S1 = S1(1 : paddedHalfLen, :);
% Optimize real matrix with extended spectral(Extended spectral simulate circular convolution)
ext = S1;
clear sg
sg.S = ext;
sg.target = target;
func3_ = @(x) corrRealMtxFastExtendedPagedHC2(x, paddedHalfLen, paddedFFTLen, prepad, pospad, hop, sg, [], synWnd2, fftLen, eachSide);
preview2_ = func3_(weightsOpt3);
function [S, t] = ltv_spectrogram3(x, coeff, cpxRes)
x = x(:);
x = [zeros(coeff.fftLen - coeff.hop, 1); x; zeros(coeff.fftLen / 2, 1)];
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.fftLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
x = [x; zeros(length(x) - paddedLen + coeff.hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
S = zeros(nframes, coeff.halfLen, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.fftLen - 1);
    X = fft(frame);
    S(i, :) = ltv_1st_2ndNoSq3(X, coeff.b, coeff.a, coeff.c1, coeff.c2, coeff.prepad, coeff.pospad, coeff.halfLen, cpxRes);
end
t = (frmIdx - 1 + coeff.fftLen / 2) / coeff.fs;
S = permute(S, [2, 1, 3]);
end
function spec = ltv_1st_2ndNoSq3(dftSpec, b, a, c1, c2, prepad, pospad, halfLen, cpxRes)
spec = cpxRes * dftSpec;
spec = spec(1 : halfLen);
spec(1, :) = real(spec(1, :));
spec(end, :) = real(spec(end, :));
end
function y = corrRealMtxFastExtendedPagedHC2(weights, paddedHalfLen, paddedFFTLen, prepad, pospad, hop, sigs, gWeights, synWnd, fftLen, eachSide) % Same as function corrRealMtxFastExtendedPaged
x = reshape(weights, paddedFFTLen, paddedFFTLen);
q_fft_frameinv2Re = pagemtimes(x, real(sigs.S));
q_fft_frameinv2Im = pagemtimes(x, imag(sigs.S));
q_fft_frameinv2Im(1, :, :) = 0;
q_fft_frameinv2Im(end, :, :) = 0;
mtx = real(ifft(q_fft_frameinv2Re + q_fft_frameinv2Im * 1j));
y3 = overlapAdd3D(mtx, hop);
y = y3(fftLen - hop + eachSide + 1 : end, :);
end