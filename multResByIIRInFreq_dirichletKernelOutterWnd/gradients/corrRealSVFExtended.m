function [f, g, outInv] = corrRealSVFExtended(weights, halfLen, fftLen, prepad, pospad, hop, sigs, gWeights)
poles = ADNode(weights);
nSections = 1;
coefNum = halfLen + prepad + pospad - 1;
b0 = reshape(poles(coefNum * nSections * 0 + 1 : coefNum * nSections * 1), [coefNum, nSections]);
d1 = reshape(poles(coefNum * nSections * 1 + 1 : coefNum * nSections * 2), [coefNum, nSections]);
d2 = reshape(poles(coefNum * nSections * 2 + 1 : coefNum * nSections * 3), [coefNum, nSections]);
c1 = reshape(poles(coefNum * nSections * 3 + 1 : coefNum * nSections * 4), [coefNum, nSections]);
c2 = reshape(poles(coefNum * nSections * 4 + 1 : coefNum * nSections * 5), [coefNum, nSections]);
for idx = 1 : size(sigs, 1)
    %% Inverse transform
    % Frequency weighted interpolation between original and windowed spectra
    q_fft_frameinv2Re = svfFilter(sigs{idx}.SRe, b0, d1, d2, c1, c2, 0, 0, 0, 0);
    q_fft_frameinv2ReRev = svfFilter(q_fft_frameinv2Re, b0, d1, d2, c1, c2, 0, 0, 0, 1);
    q_fft_frameinv2Im = svfFilter(sigs{idx}.SIm, b0, d1, d2, c1, c2, 0, 0, 0, 0);
    q_fft_frameinv2ImRev = svfFilter(q_fft_frameinv2Im, b0, d1, d2, c1, c2, 0, 0, 0, 1);
    q_fft_frameinv2ReCut = q_fft_frameinv2ReRev(prepad + 1 : end - pospad + 1, :);
    q_fft_frameinv2ImCut = q_fft_frameinv2ImRev(prepad + 1 : end - pospad + 1, :);
    td = separate_irfft(q_fft_frameinv2ReCut, q_fft_frameinv2ImCut, fftLen);
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
if ~isempty(gWeights)
    g = g .* gWeights;
end
g = g(:);
end