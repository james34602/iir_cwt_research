rng(1)
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
phaseShifter4 = 0.5 * exp(-1i * (2*pi*(0:halfLen-1)/fftLen).');
phaseShifter5 = 0.5 * exp(1i * (2*pi/fftLen).');
phaseShifter6 = 0.5 * exp(-1i * (2*pi/fftLen).');
fWnd2 = phaseShifter1 .* (dftSpec1 - phaseShifter2 .* [conj(dftSpec1(2)); dftSpec1(1 : end - 1)] - phaseShifter3 .* [dftSpec1(2 : end); conj(dftSpec1(end - 1))]);
fWnd2(end) = real(fWnd2(end));
sig3 = ifft(createSym(fWnd2));
fWnd3 = phaseShifter4 .* (dftSpec1 - phaseShifter5 .* [conj(dftSpec1(2)); dftSpec1(1 : end - 1)] - phaseShifter6 .* [dftSpec1(2 : end); conj(dftSpec1(end - 1))]);
fWnd3(end) = real(fWnd3(end));
sig4 = ifft(createSym(fWnd3));
figure(1)
plot(sig2);
hold on
plot(sig3);
plot(sig4);
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
fWnd5 = phaseShifter4 .* (dftSpec2 + phaseShifter5 .* [conj(dftSpec2(2)); dftSpec2(1 : end - 1)] + phaseShifter6 .* [dftSpec2(2 : end); conj(dftSpec2(end - 1))]);
fWnd5(end) = real(fWnd5(end));
sig5 = ifft(createSym(fWnd4));
sig6 = ifft(createSym(fWnd5));
figure(2)
plot(ifftshift(sig4))
hold on
plot(ifftshift(sig5))
plot(ifftshift(sig6))
hold off
axis tight;
function y = createSym(x)
halfLen = length(x);
fftLen = (halfLen - 1) * 2;
y = x;
y(halfLen+1:fftLen,:) = conj(y(halfLen-1:-1:2,:));
end