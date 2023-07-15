function [f, g, outInv] = corrRealMtxFastExtended(weights, halfLen, fftLen, prepad, pospad, hop, sigs, gWeights)
x = ADNode(reshape(weights, halfLen, halfLen + prepad + pospad - 1));
for idx = 1 : size(sigs, 1)
    %% Inverse transform
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2Re = x * sigs{idx}.SRe;
    q_fft_frameinv2Im = x * sigs{idx}.SIm;
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
f = mse.value;
g = mse.backprop(1);
g = g(:);
if ~isempty(gWeights)
    g = g .* gWeights;
end
end