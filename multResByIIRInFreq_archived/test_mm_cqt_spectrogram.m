nfft = 4096;
hop = 512;
[y, fs] = audioread('../sinesweep.wav');

disp('Calculating STFT Spectrogram.')
[s,f,t,p] = spectrogram([y; zeros(nfft, 1)],hann(nfft),nfft-hop,nfft,fs); 

Noct = 28;
disp('Calculating IIR CQT Spectrogram.')
tic
[s_mmq, f_mmq, t_mmq] = mm_cqt_complex(y,nfft,hop,fs,Noct);
toc

%% plots

bias = 0.001;

stft_dB = 20*log10(abs(s)+bias);
mmcqt_stft_dB = 20*log10(abs(s_mmq)+bias);

max_freq = 20000;

[m,kmax_sfft] = min(abs(f-max_freq));
[m,kmax_mmcqt] = min(abs(f_mmq-max_freq));

figure(1)
dim = size(stft_dB);
imagesc(t,f(1:kmax_sfft),-stft_dB(1:kmax_sfft,:));
colormap(jet);
set(gca,'YDir','normal')
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('STFT Spectrogram (magnitude in dB)')

figure(2)
imagesc(t_mmq,f_mmq(1:kmax_mmcqt),-mmcqt_stft_dB(1:kmax_mmcqt,:))
colormap(jet);
set(gca,'YDir','normal')
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('MM Constant Q Spectrogram (magnitude in dB)')