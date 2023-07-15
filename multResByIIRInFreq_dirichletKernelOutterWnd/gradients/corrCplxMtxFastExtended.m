function [f, g, outInv] = corrCplxMtxFastExtended(weights, halfLen, fftLen, prepad, pospad, hop, sigs, gWeights)
x = ADNode(reshape(weights, (halfLen + prepad + pospad - 1), (halfLen + prepad + pospad - 1) * 2));
xRe = x(:, 1 : (halfLen + prepad + pospad - 1));
xIm = x(:, (halfLen + prepad + pospad - 1) + 1 : (halfLen + prepad + pospad - 1) * 2);
for idx = 1 : size(sigs, 1)
    q_fft_frameinv2Re = xRe * sigs{idx}.SRe - xIm * sigs{idx}.SIm;
    q_fft_frameinv2Im = xRe * sigs{idx}.SIm + xIm * sigs{idx}.SRe;
    q_fft_frameinv2ReCut = q_fft_frameinv2Re(prepad + 1 : end - pospad + 1, :);
    q_fft_frameinv2ImCut = q_fft_frameinv2Im(prepad + 1 : end - pospad + 1, :);
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
g = g(:);
if ~isempty(gWeights)
    g = g .* gWeights;
end
end