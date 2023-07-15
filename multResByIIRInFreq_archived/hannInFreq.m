rng(1)
fftLen = 1024;
actualWnd = hann(fftLen, 'periodic');
halfLen = fftLen / 2 + 1;
sig = (rand(fftLen, 1) - 0.5) * 2;
% sig = ones(fftLen, 1);
%% Hann in frequency domain without input being shifted by fftshift
x_fft1 = fft(sig);
fWnd1 = 0.5 * (x_fft1 - 0.5 * circshift(x_fft1, 1) - 0.5 * circshift(x_fft1, -1));
sig2 = real(ifft(fWnd1));
%% Hann in frequency domain without 1 sample delayed input being shifted by fftshift
phaseShifter1 = 0.5 * exp(1i * (2*pi*(0:fftLen-1)/fftLen).');
phaseShifter2 = 0.5 * exp(-1i * (2*pi/fftLen).');
phaseShifter3 = 0.5 * exp(1i * (2*pi/fftLen).');
phaseShifter4 = 0.5 * exp(-1i * (2*pi*(0:fftLen-1)/fftLen).');
phaseShifter5 = 0.5 * exp(1i * (2*pi/fftLen).');
phaseShifter6 = 0.5 * exp(-1i * (2*pi/fftLen).');
fWnd2 = phaseShifter1 .* (x_fft1 - phaseShifter2 .* circshift(x_fft1, 1) - phaseShifter3 .* circshift(x_fft1, -1));
sig3 = real(ifft(fWnd2));
fWnd3 = phaseShifter4 .* (x_fft1 - phaseShifter5 .* circshift(x_fft1, 1) - phaseShifter6 .* circshift(x_fft1, -1));
sig4 = real(ifft(fWnd3));
figure(1)
plot(sig2);
hold on
plot(sig3);
plot(sig4);
hold off
axis tight;
%% Prepare a fftshifted input
x_fft2 = fft(fftshift(sig));
%% Hann in frequency domain with input being shifted by fftshift
q_fft_frameinv = 0.5 * (x_fft2 + 0.5 * circshift(x_fft2, 1) + 0.5 * circshift(x_fft2, -1));
sig4 = real(ifft(q_fft_frameinv));
%% Hann in frequency domain with 1 sample delayed input being shifted by fftshift
q_fft_frameinv2 = phaseShifter1 .* (x_fft2 + phaseShifter2 .* circshift(x_fft2, 1) + phaseShifter3 .* circshift(x_fft2, -1));
sig5 = real(ifft(q_fft_frameinv2));
q_fft_frameinv3 = phaseShifter4 .* (x_fft2 + phaseShifter5 .* circshift(x_fft2, 1) + phaseShifter6 .* circshift(x_fft2, -1));
sig6 = real(ifft(q_fft_frameinv3));
figure(2)
plot(ifftshift(sig4))
hold on
plot(ifftshift(sig5))
plot(ifftshift(sig6))
hold off
axis tight;