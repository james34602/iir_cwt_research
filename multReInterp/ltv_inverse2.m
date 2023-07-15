function y = ltv_inverse2(S,fftLen,hop,reqSynthesisWnd,correctionWndHF,wndCorrectionWeightingLF,wndCorrectionWeightingHF)
nframes = size(S, 2);
halfLen = fftLen / 2 + 1;
accFrame2 = zeros(fftLen, 1);
frameindex = 1 + (0 : nframes - 1) * hop;
y = zeros(nframes * hop, 1);
for i=1:nframes
    %% Inverse transform
    q_fft_frameinv = S(:, i);
    if ~isempty(wndCorrectionWeightingLF)
        q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
        correctedTime = ifft(q_fft_frameinv);
        correctedTime2 = correctedTime .* correctionWndHF;
        getbackCorrectedToSpectrum2 = fft(correctedTime2);
        % Frequency weighted interpolation between original and windowed spectra
        q_fft_frameinv = S(:, i) .* wndCorrectionWeightingLF + getbackCorrectedToSpectrum2(1 : halfLen) .* wndCorrectionWeightingHF;
    end
    if reqSynthesisWnd
        % Virtually multiplying Hann window in time domain on frequency domain
        q_fft_frameinv = 0.25 * (2 * q_fft_frameinv + [conj(q_fft_frameinv(2)); q_fft_frameinv(1 : end - 1)] + [q_fft_frameinv(2 : end); conj(q_fft_frameinv(end - 1))]);
    end
    % Reflect spectrum and conjugate
    q_fft_frameinv(halfLen+1:fftLen,:) = conj(q_fft_frameinv(halfLen-1:-1:2,:));
    yInvShifted = ifftshift(ifft(q_fft_frameinv));
    % Overlap add
    accFrame2 = accFrame2 + yInvShifted;
    myOut = accFrame2(1 : hop);
    accFrame2 = [accFrame2(hop + 1 : end); zeros(hop, 1)];
    y(frameindex(i):frameindex(i)+hop-1) = myOut;
end
end