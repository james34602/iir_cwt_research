% rng(1)
fftLen = 1024;
actualWnd = hann(fftLen, 'periodic');
halfLen = fftLen / 2 + 1;
sig = (rand(fftLen, 1) - 0.5) * 2;
% sig = ones(fftLen, 1);
%% Hann in frequency domain without input being shifted by fftshift
x_fft1 = fft(sig);
dftSpec1 = x_fft1(1 : halfLen);
fWnd1 = 0.5 * (dftSpec1 - 0.5 * [conj(dftSpec1(2)); dftSpec1(1 : end - 1)] - 0.5 * [dftSpec1(2 : end); conj(dftSpec1(end - 1))]);
sig2 = ifft(createSym(fWnd1));
%% Hann in frequency domain without 1 sample delayed input being shifted by fftshift
phaseShifter1 = 0.5 * exp(1i * (2*pi*(0:halfLen-1)/fftLen).');
phaseShifter2 = 0.5 * exp(-1i * (2*pi/fftLen).');
phaseShifter3 = 0.5 * exp(1i * (2*pi/fftLen).');
fWnd2 = phaseShifter1 .* (dftSpec1 - phaseShifter2 .* [conj(dftSpec1(2)); dftSpec1(1 : end - 1)] - phaseShifter3 .* [dftSpec1(2 : end); conj(dftSpec1(end - 1))]);
fWnd2(end) = real(fWnd2(end));
sig3 = ifft(createSym(fWnd2));
figure(1)
plot(sig2 - sig .* actualWnd);
hold on
plot(sig3 - circshift(sig, -1) .* actualWnd);
hold off
axis tight;
%% Prepare a fftshifted input
x_fft2 = fft(fftshift(sig));
dftSpec2 = x_fft2(1 : halfLen);
%% Hann in frequency domain with input being shifted by fftshift
fWnd3 = 0.5 * (dftSpec2 + 0.5 * [conj(dftSpec2(2)); dftSpec2(1 : end - 1)] + 0.5 * [dftSpec2(2 : end); conj(dftSpec2(end - 1))]);
sig4 = ifft(createSym(fWnd3));
%% Hann in frequency domain with 1 sample delayed input being shifted by fftshift
fWnd4 = phaseShifter1 .* (dftSpec2 + phaseShifter2 .* [conj(dftSpec2(2)); dftSpec2(1 : end - 1)] + phaseShifter3 .* [dftSpec2(2 : end); conj(dftSpec2(end - 1))]);
fWnd4(end) = real(fWnd4(end));
sig5 = ifft(createSym(fWnd4));
figure(2)
plot(ifftshift(sig4) - sig .* actualWnd)
hold on
plot(ifftshift(sig5) - circshift(sig, -1) .* actualWnd)
hold off
axis tight;
function y = createSym(x)
halfLen = length(x);
fftLen = (halfLen - 1) * 2;
y = x;
y(halfLen+1:fftLen,:) = conj(y(halfLen-1:-1:2,:));
end