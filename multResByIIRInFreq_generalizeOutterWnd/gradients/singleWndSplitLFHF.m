function [f, g] = singleWndSplitLFHF(weights, halfLen, frameSize, hop, sigs)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingLF(halfLen+1:frameSize) = wndCorrectionWeightingLF(halfLen-1:-1:2);
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
wndCorrectionWeightingHF(halfLen+1:frameSize) = wndCorrectionWeightingHF(halfLen-1:-1:2);
for idx = 1 : size(sigs, 1)
    %% Inverse transform
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2 = sigs{idx}.S .* wndCorrectionWeightingLF + sigs{idx}.getbackCorrectedToSpectrum .* wndCorrectionWeightingHF;
    mtx = circshift(ifft(q_fft_frameinv2), frameSize / 2);
    y3 = fold(mtx, frameSize, hop);
    y = y3(frameSize - hop + 1 : frameSize - hop + frameSize);
    meaSqr = (y - sigs{idx}.target) .^ 2;
    if idx == 1
        mse = mean(meaSqr);
    else
        mse = mse + mean(meaSqr);
    end
end
mse = mse / size(sigs, 1);
if nargout > 0
    f = mse.value;
end
if nargout > 1
    g = real(mse.backprop(1));
end
end