function [f, g, outInv] = singleWndSplitLFHFMtx(weights, fftLen, prepad, pospad, halfLen, hop, sigs, weightsOpt)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
extLF = vertcat(vertcat(wndCorrectionWeightingLF(prepad + 1 : -1 : 2, :), wndCorrectionWeightingLF), wndCorrectionWeightingLF(halfLen - 1 : -1 : (halfLen - pospad + 1), :));
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
extHF = vertcat(vertcat(wndCorrectionWeightingHF(prepad + 1 : -1 : 2, :), wndCorrectionWeightingHF), wndCorrectionWeightingHF(halfLen - 1 : -1 : (halfLen - pospad + 1), :));
for idx = 1 : size(sigs, 1)
    %% Inverse transform
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2Re = sigs{idx}.SRe .* extLF + sigs{idx}.getbackCorrectedToSpectrumRe .* extHF;
    q_fft_frameinv2Im = sigs{idx}.SIm .* extLF + sigs{idx}.getbackCorrectedToSpectrumIm .* extHF;
    q_fft_frameinv2Re = weightsOpt * q_fft_frameinv2Re;
    q_fft_frameinv2Im = weightsOpt * q_fft_frameinv2Im;
    q_fft_frameinv2Im(1, :) = 0;
    q_fft_frameinv2Im(end, :) = 0;
    td = separate_irfft(q_fft_frameinv2Re, q_fft_frameinv2Im, fftLen);
    mtx = circshift(td, fftLen / 2);
    y3 = fold(mtx, fftLen, hop);
    y = y3(fftLen - hop + 1 : size(sigs{idx}.target, 1) - hop + fftLen);
    meaSqr = (y - sigs{idx}.target) .^ 2;
    if idx == 1
        mse = mean(meaSqr);
        outInv = y.value;
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