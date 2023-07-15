function [S, t] = iirltv_spectrogram(x, coeff)
x = [x(:); zeros(coeff.halfLen - 2, 1)];
ny = length(x);
% zero padding at the end to complete the last frame
frmIdx = 1 + (0 : ny - 1) * coeff.hop;
% matrix to store the complex spectrogram
S = (zeros(ny, coeff.fftLen/2+1, 2, 'like', 1i));
shiftCentre = zeros(coeff.halfLen, 1);
for i = 1 : coeff.halfLen
    if ~bitget(i, 1)
        shiftCentre(i) = -1;
    else
        shiftCentre(i) = 1;
    end
end
W = exp((0 : (coeff.halfLen - 1))' * (2 * pi) / coeff.fftLen * 1i);
fifo = zeros(coeff.fftLen, 1);
sdft = zeros(coeff.halfLen, 1);
xn_M = 0;
%% IIR LTV Q FFT transform
for i=1:ny
    fifo = [fifo(2 : coeff.fftLen); x(i)];
    sdft = W .* (sdft - xn_M + x(i));
    %% new xn_M
    xn_M = fifo(1);
    S(i, :, :) = ltv_1st_2nd(sdft .* shiftCentre, coeff.b, coeff.a, coeff.c1, coeff.c2, coeff.prepad, coeff.pospad, coeff.phaseShifter1, coeff.phaseShifter2, coeff.phaseShifter3, ...
        coeff.corrF, coeff.fftLen, coeff.halfLen, 1);
end
t = (frmIdx - 1 + coeff.fftLen / 2) / coeff.fs;
S = permute(S, [2, 1, 3]);
end