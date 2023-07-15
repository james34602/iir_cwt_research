function [f, g, outInv] = corrCplxMtx(weights, halfLen, fftLen, hop, sigs, WinvdftRe, WinvdftIm, gWeights)
x = ADNode(reshape(weights, halfLen, halfLen * 2));
xRe = x(:, 1 : halfLen);
xIm = x(:, halfLen + 1 : halfLen * 2);
for idx = 1 : size(sigs, 1)
    q_fft_frameinv2Re = xRe * sigs{idx}.SRe - xIm * sigs{idx}.SIm;
    q_fft_frameinv2Im = xRe * sigs{idx}.SIm + xIm * sigs{idx}.SRe;
    if mod(fftLen, 2) == 0
        q_fft_frameinv2Re(halfLen+1:fftLen,:) = q_fft_frameinv2Re(halfLen-1:-1:2,:);
        q_fft_frameinv2Im(halfLen+1:fftLen,:) = q_fft_frameinv2Im(halfLen-1:-1:2,:);
    else
        q_fft_frameinv2Re(halfLen+1:fftLen,:) = q_fft_frameinv2Re(halfLen:-1:2,:);
        q_fft_frameinv2Im(halfLen+1:fftLen,:) = q_fft_frameinv2Im(halfLen:-1:2,:);
    end
    td = WinvdftRe * q_fft_frameinv2Re - WinvdftIm * q_fft_frameinv2Im;
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
% g = [g(:, 1 : halfLen)', g(:, halfLen + 1 : halfLen * 2)'];
g = g(:);
end