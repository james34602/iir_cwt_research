function [f, g, outInv] = corrRealMtxFastExtendedPaged(weights, halfLen, fftLen, prepad, pospad, hop, sigs, gWeights)
x = ADNode(reshape(weights, halfLen, halfLen + prepad + pospad - 1));
q_fft_frameinv2Re = pagemtimes(x, sigs.SRe);
q_fft_frameinv2Im = pagemtimes(x, sigs.SIm);
q_fft_frameinv2Im(1, :, :) = 0;
q_fft_frameinv2Im(end, :, :) = 0;
td = separate_irfft(q_fft_frameinv2Re, q_fft_frameinv2Im, fftLen);
mtx = circshiftCustom(td, fftLen / 2);
y3 = fold3D(mtx, fftLen, hop);
y = y3(fftLen - hop + 1 : size(sigs.target, 1) - hop + fftLen, :);
outInv = y.value(:, 1);
mse = (y - sigs.target) .^ 2;
mse = mean(mse(:));
f = mse.value;
g = mse.backprop(1);
g = g(:);
if ~isempty(gWeights)
    g = g .* gWeights;
end
end
function b = circshiftCustom(a,p)
m = size(a,1);
b = a(mod((0:m-1)-double(rem(p,m)), m)+1, :, :);
end