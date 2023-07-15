frameSize = 2048;
halfLen = frameSize / 2 + 1;
hop = 32;
Q = 3;
% numberOfVirtualBands = 10 * oct + 1
% Q = (sqrt(2 ^ (1 / oct)) / ((2 ^ (1 / oct)) - 1))
% disp('Calculating STFT Spectrogram.')
% [s,f,t,p] = spectrogram([y; zeros(nfft, 1)],hann(nfft),nfft-hop,nfft,fs);
[y, fs] = loadSignal(7, frameSize);
disp('Calculating IIR CQT Spectrogram.')
tic
[s_q,f_q,t_q] = iir_cqt_spectrogram(y,frameSize,hop,fs,Q);
toc
% Normalize to 0dB
s_q = s_q ./ max(abs(s_q(:)));
max(20 * log10(abs(s_q)), [], 'all')
% %% Inspect undersampling phenomenon by plotting sine sweep spectrum slice
% normalizePeak = max(abs(s_q), [], 1);
% plot(abs(s_q) ./ normalizePeak);
% %% Inspect undersampling phenomenon by plotting sum of sine sweep spectrum slice in time
% pkNormalizedMagnitude = abs(s_q) ./ normalizePeak;
% plot(sum(pkNormalizedMagnitude, 1));
% %% Plot sum of spectrum slice in time
% pkNormalizedMagnitude = abs(s_q);
% plot(sum(pkNormalizedMagnitude, 1));
% plot(20 * log10(abs(s_q)))
%% plots
% stft_dB = 20*log10(abs(s));
cqt_stft_dB = 20*log10(abs(s_q));

% figure(1)
% dim = size(stft_dB);
% imagesc(t,f,-stft_dB);
% colormap(jet);
% set(gca,'YDir','normal')
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('STFT Spectrogram (magnitude in dB)')

cqt_stft_dB=cqt_stft_dB(:, 1 : size(t_q, 2));
figure(2)
imagesc(t_q,f_q,cqt_stft_dB)
colormap(jet);
caxis([-100, 0])
set(gca,'YDir','normal')
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('IIR Constant Q Spectrogram (magnitude in dB)')
colorbar