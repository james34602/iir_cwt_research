function [S, t] = ltv_spectrogram(x, coeff)
x = x(:);
x = [zeros(coeff.fftLen - coeff.hop, 1); x; zeros(coeff.fftLen / 2, 1)];
ny = length(x);
nframes = fix(ceil(ny/coeff.hop) - coeff.fftLen / coeff.hop / 2);
% zero padding at the end to complete the last frame
paddedLen = coeff.hop * nframes;
x = [x; zeros(length(x) - paddedLen + coeff.hop * 2, 1)];
frmIdx = 1 + (0 : nframes - 1) * coeff.hop;
% matrix to store the complex spectrogram
S = zeros(nframes, coeff.fftLen/2+1, 3, 'like', 1i);
%% IIR LTV Q FFT transform
for i=1:nframes
    frame = x(frmIdx(i) : frmIdx(i) + coeff.fftLen - 1);
    X = fft(fftshift(frame .* coeff.analysisWnd));
    S(i, :, :) = ltv_1st_2nd(X(1 : coeff.halfLen), coeff.b, coeff.a, coeff.c1, coeff.c2, coeff.prepad, coeff.pospad, coeff.phaseShifter1, ...
        coeff.corrF, coeff.fftLen, coeff.halfLen, 1);
end
t = (frmIdx - 1 + coeff.fftLen / 2) / coeff.fs;
S = permute(S, [2, 1, 3]);
end