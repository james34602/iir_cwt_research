function [f, wt_grad, y] = corrRealMtxFastExtendedPagedHC(weights, halfLen, fftLen, prepad, pospad, hop, sigs, gWeights) % Same as function corrRealMtxFastExtendedPaged
x = reshape(weights, halfLen, halfLen + prepad + pospad - 1);
q_fft_frameinv2Re = pagemtimes(x, sigs.SRe);
q_fft_frameinv2Im = pagemtimes(x, sigs.SIm);
q_fft_frameinv2Im(1, :, :) = 0;
q_fft_frameinv2Im(end, :, :) = 0;
if mod(fftLen,2) == 0
    fg = 0;
    q_fft_frameinv2Re(halfLen+1:fftLen, :, :) = q_fft_frameinv2Re(halfLen-1:-1:2, :, :);
    q_fft_frameinv2Im(halfLen+1:fftLen, :, :) = -q_fft_frameinv2Im(halfLen-1:-1:2, :, :);
else
    fg = 1;
    q_fft_frameinv2Re(halfLen+1:fftLen, :, :) = q_fft_frameinv2Re(halfLen:-1:2, :, :);
    q_fft_frameinv2Im(halfLen+1:fftLen, :, :) = -q_fft_frameinv2Im(halfLen:-1:2, :, :);
end
mtx = circshiftCustom(ifft(q_fft_frameinv2Re + q_fft_frameinv2Im * 1j), fftLen / 2);
y3 = overlapAdd3D(mtx, hop);
y = y3(fftLen - hop + 1 : size(sigs.target, 1) - hop + fftLen, :);
dif = y - sigs.target;
avgG = 1 / (size(sigs.target, 1) * size(sigs.target, 2));
f = sum(dif(:) .^ 2) * avgG;
se_grad = ones(size(sigs.target)) * avgG .* dif * 2;
y_grad = zeros(size(y3));
y_grad(fftLen - hop + 1 : size(sigs.target, 1) - hop + fftLen, :) = se_grad;
rec = zeros(fftLen, ceil(size(y_grad, 1) / hop), size(y_grad, 2));
for idx = 1 : size(y_grad, 2)
    rec(:, :, idx) = buffer(y_grad(:, idx), fftLen, fftLen - hop);
end
tmp = ifft(circshiftCustom(rec(:, fftLen / hop : end, :), fftLen / 2));
if fg == 0
    tmp = tmp(1 : halfLen, :, :);
    tmp(2 : end - 1, :, :) = tmp(2 : end - 1, :, :) * 2;
else
    tmp = tmp(1 : halfLen, :, :);
    tmp(2 : end, :, :) = tmp(2 : end, :, :) * 2;
end
xRe_grad = real(tmp);
xIm_grad = -imag(tmp);
xIm_grad(1, :, :) = 0;
xIm_grad(end, :, :) = 0;
wt_grad = sum(pagemtimes(xIm_grad, 'none', sigs.SIm, 'transpose'), 3) + sum(pagemtimes(xRe_grad, 'none', sigs.SRe, 'transpose'), 3);
wt_grad = wt_grad(:);
if ~isempty(gWeights)
    wt_grad = wt_grad .* gWeights;
end
end
function b = circshiftCustom(a,p)
m = size(a,1);
b = a(mod((0:m-1)-double(rem(p,m)), m)+1, :, :);
end