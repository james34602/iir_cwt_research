rng(1)
fftLen = 16;
actualWnd = hann(fftLen, 'periodic');
halfLen = fftLen / 2 + 1;
% sig = (rand(fftLen, 1) - 0.5) * 2;
sig = ones(fftLen, 1);
%% Hann in frequency domain without input being shifted by fftshift
dftM = dftmtx(fftLen);
freq1 = sig .* conj(dftM(:, 5)) * 0.5;
freq2 = sig .* conj(dftM(:, 6)) * 0.5;
x_fft1 = fft(freq1 + freq2);
freq1 = hann(fftLen, 'periodic') .* freq1;
freq2 = hann(fftLen, 'periodic').^2 .* freq2;
x_fft2 = fft(freq1 + freq2);x_fft2(abs(x_fft2) < eps*20) = 0;x_fft2 = real(x_fft2);
x_fft1(abs(x_fft1) < eps*20) = 0;x_fft1 = real(x_fft1);
% [seg1, rng1] = wndingAutoPad(x_fft1, 4 : 6);
% [seg2, rng2] = wndingAutoPad(x_fft1, 5 : 7);
% [~, rng3] = wndingAutoPad(x_fft1, rng2);
% [seg3, ~] = wndingAutoPad(seg2, 1 : length(seg2));
% x_fft3(rng1) = x_fft3(rng1) + seg1;
%  x_fft3(rng2) = x_fft3(rng2) + seg2;
% % x_fft3(rng3) = x_fft3(rng3) + seg3;
kernel = [-1, 2, -1] / 4;
% ker = circshift([ker; zeros(floor((fftLen - 3)/1), 1)], -1);
% ker = circshift(ker, -1);
ker = zeros(fftLen, 3);
ker(:, 2) = 1;
ker(5 : 6, :) = [kernel; kernel];
% ker(6, 1 : 3) = [0, 1, 0];
% x_fft3 = circular_conv(ker, x_fft1);
x_fft3 = dconv(x_fft1, ker);
x_fft3 = x_fft3(2 : end - 1);
ker(7, 1 : 3) = kernel;
x_fft3 = dconv(x_fft3, ker);
x_fft3 = x_fft3(2 : end - 1);
fWnd1 = wnding(x_fft1);
fWnd1 = [fWnd1, wnding(fWnd1)];
sig2 = ifft(fWnd1);
% figure(1)
plot(real(sig2));
hold on
plot(imag(sig2));
hold off
axis tight;
function vc = dconv(vs, h)
N = length(vs);
M = size(h, 2);   
lout=N+M-1;
vc=zeros(1,lout); 
for i = 1:N
    vc(i + (1 : M) - 1) = vc(i + (1 : M) - 1) + h(i, :) .* vs(i);
end
end