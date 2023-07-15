frameSize = 512;
halfLen = frameSize / 2 + 1;
hop = 1;
oct = 5;
HFSamplingLimit = 0.7;
order = 2;
reqSynthesisWnd = 1;
% numberOfVirtualBands = 10 * oct + 1
% Q = (sqrt(2 ^ (1 / oct)) / ((2 ^ (1 / oct)) - 1))
% disp('Calculating STFT Spectrogram.')
[x, fs] = loadSignal(6, frameSize);
% x = x(1 : 4096);
% kd = zeros(80, 1);
% kd(10) = 15;
% kd(50) = 15;
% 
% fs = 2; % Sampling frequency (samples per second)
% dt = 1/fs; % seconds per sample
% StopTime = 210; % Samples
% t = (0:dt:StopTime)'; % seconds
% F2 = 0.007; % Sine wave frequency (hertz)
% depth = 5;
% 
% F = 0.75; % Sine wave frequency (hertz)
% y = sin(2*pi*F.* (t + (1 + cos(2*pi.*F2.*t + pi / 4) * depth) - 1));
% comp1 = gradient(F.* (t + (1 + cos(2*pi.*F2.*t + pi / 4) * depth) - 1));
% comp1 = comp1 - mean(comp1);
% plot(comp1)
% 
% D = length(y) / fs; % duration in seconds
% f1 = 0.25;
% f2 = 0.6;
% bandlimitChirp = chirp2(t, D, f1, f2);
% y = y + bandlimitChirp;
% F = 0.2; % Sine wave frequency (hertz)
% x = [kd; y + sin(2*pi*F.* t)];
% 
% 
% N  = 1024;
% t  = (0:N-1)/sqrt(N);
% %% Test signal 1
% s1 = sin(2*pi*(250*t/sqrt(N)+50*(t/sqrt(N)).^3));
% s2 = sin(2*pi*(130*t/sqrt(N)+100*(t/sqrt(N)).^2));
% s3 = sin(2*pi*(90*t/sqrt(N)+0.2*cos(3*pi*(t/sqrt(N)))));
% s = s1+s2+s3;
% %% test signal 2
% s1 = sin(2*pi*(330*t/sqrt(N)+16*cos(3*pi*t/sqrt(N))));
% s2 = sin(2*pi*(190*t/sqrt(N)+9*cos(3*pi*t/sqrt(N))));
% s3 = sin(2*pi*(40*t/sqrt(N)));
% s_ = s1+s2+s3;
% 
% x = s';
[coeff, f_q] = ltv_precomute2(frameSize, hop, fs, oct, order, HFSamplingLimit, 1, 1, reqSynthesisWnd);
disp('Calculating IIR CQT Spectrogram.')
e = norm(x,'fro').^2;
if any(isnan(x)) || any(isinf(x))
    disp('NaN or Inf in signal is not allowed, reassignment algorithm will crash');
end
t_q = 0 : (length(x) - 1);
tic
spec = ltv_spectrogram(x, coeff);
toc
figure(1)
[matSpec, cwtF] = wsst(x, fs, 'VoicesPerOctave', 48, 'amor');
s_q = matSpec ./ max(abs(matSpec(:)));
h = pcolor(t_q, cwtF, 20 * log10(abs(s_q)));
set(h,'EdgeColor', 'none')
colorbar
caxis([-100, 0])
title('WSST')
% figure(2)
% [imf,residual,info] = vmd(x);
% matSpec=hht(imf);
% s_q = matSpec ./ max(abs(matSpec(:)));
% cwtF = linspace(0, 1, size(matSpec, 1));
% h = pcolor(0 : (length(x) - 1), cwtF, 20 * log10(abs(s_q)));
% set(h,'EdgeColor', 'none')
% colorbar
% caxis([-100, 0])
% max(20 * log10( abs(spec(:, :, 1)) * 2 ./ (sum(hann(frameSize, 'periodic'))) ), [], 'all')
inverseSynSq = 1;
ridgeExt = 0;
spec2 = spec(:, :, 1);
spec2(halfLen+1:frameSize,:) = conj(spec2(halfLen-1:-1:2,:));
E3 = norm(abs(spec2),'fro').^2;
synSq = spec(:, frameSize / 2 + 0 : end, 2);
synSq = synSq(:, 1 : length(x));
% Normalize to 0dB
s_q = synSq ./ max(abs(synSq(:)));
%% Inverse SST
if inverseSynSq && hop == 1
    if ridgeExt
        contour(t_q, f_q, abs(synSq))
        [fridge,iridge] = wsstridge(synSq,5,f_q(:),'NumRidges',3, 'NumFrequencyBins', 24);
        hold on;
        plot(t_q,fridge,'k--','linewidth',2);
        hold off;

        nbins = 4;
        sizIridge = size(iridge);
        xrec = zeros(sizIridge,'like',synSq);
        for j = 1:sizIridge(1,2)
            Mask = zeros(size(synSq),'like',synSq);
            for i = 1:sizIridge(1,1)
                Mask(max(iridge(i,j)-nbins(1),1):min(iridge(i,j)+nbins(1),size(synSq,1)),i) = 1;
            end
            % Do the inverse
            maskSST = synSq.*Mask;
            xrec(:,j) = real(sum(maskSST,1) + conj(sum(maskSST(2 : frameSize / 2,:), 1)))';
        end
    else
        xrec = real(sum(synSq, 1) + conj( sum(synSq(2 : frameSize / 2, :), 1) ));
        %xrec = iwsst(matSpec);
        SNR = 10*log10(sum(abs(x).^2)/sum(abs(x-xrec(:)).^2)) %calculate the SNR
    end
    figure(3)
    plot(xrec);
    axis tight;
end
%
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
% axis tight
% ylim([-30, 0])

%% Plot spectrogram
figure(4)
imagesc(t_q, f_q, 20 * log10(abs(s_q)))
caxis([-100, 0])
colormap(jet);
set(gca,'YDir','normal');
colorbar
if order == 2
    title('My second order')
else
    title('My first order')
end
function x = chirp2(t,t1,f0,f1)
beta = (f1-f0)./t1;
x = cos(2*pi * ( 0.5* beta .* (t .* t) + f0 * t));
end