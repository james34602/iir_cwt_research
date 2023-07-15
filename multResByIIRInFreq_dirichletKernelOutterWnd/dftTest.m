rng(1)
fftLen = 77;
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
Wdft = fft(eye(fftLen));
WdftRe = real(Wdft(1 : halfLen, :));
WdftIm = imag(Wdft(1 : halfLen, :));
Winvdft = conj(Wdft) / fftLen;
WinvdftRe = real(Winvdft);
WinvdftIm = imag(Winvdft); WinvdftIm(:, halfLen + 1 : end) = -WinvdftIm(:, halfLen + 1: end);
clear Wdft Winvdft

sig = rand(fftLen, 6);
standardFFT = fft(sig);
standardFFT = standardFFT(1 : halfLen, :);

specRe = WdftRe * sig;
specIm = WdftIm * sig;
disp(mean(abs( (specRe + specIm * 1j) - standardFFT)))
if mod(fftLen, 2) == 0
    specRe(halfLen+1:fftLen,:) = specRe(halfLen-1:-1:2,:);
    specIm(halfLen+1:fftLen,:) = specIm(halfLen-1:-1:2,:);
else
    specRe(halfLen+1:fftLen,:) = specRe(halfLen:-1:2,:);
    specIm(halfLen+1:fftLen,:) = specIm(halfLen:-1:2,:);
end
td = WinvdftRe * specRe - WinvdftIm * specIm;
disp(mean(abs(td - sig)))