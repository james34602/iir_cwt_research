function [f, g] = singleWndSplitLFHFPaged(weights, halfLen, fftLen, hop, sigs)
x = ADNode(weights);
wndCorrectionWeightingLF = x(1 : halfLen);
wndCorrectionWeightingHF = x(halfLen + 1 : halfLen + halfLen);
q_fft_frameinv2Re = sigs.SRe .* wndCorrectionWeightingLF + sigs.getbackCorrectedToSpectrumRe .* wndCorrectionWeightingHF;
q_fft_frameinv2Im = sigs.SIm .* wndCorrectionWeightingLF + sigs.getbackCorrectedToSpectrumIm .* wndCorrectionWeightingHF;
td = separate_irfft(q_fft_frameinv2Re, q_fft_frameinv2Im, fftLen);
mtx = circshiftCustom(td, fftLen / 2);
y3 = fold3D(mtx, fftLen, hop);
y = y3(fftLen - hop + 1 : size(sigs.target, 1) - hop + fftLen, :);
mse = (y - sigs.target) .^ 2;
mse = mean(mse(:));
f = mse.value;
g = mse.backprop(1);
g = g(:);
end
function b = circshiftCustom(a,p)
m = size(a,1);
b = a(mod((0:m-1)-double(rem(p,m)), m)+1, :, :);
end