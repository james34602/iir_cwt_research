dim1 = 3; dim2 = 4; %[5, 2], [6, 3]
fftLen = coeffCollection{dim1, dim2}.fftLen;
hop = coeffCollection{dim1, dim2}.hop;
oct = coeffCollection{dim1, dim2}.oct;
halfLen = fftLen / 2 + 1;
HFSamplingLimit = coeffCollection{dim1, dim2}.HFSamplingLimit;
order = 2;
% numberOfVirtualBands = 10 * oct + 1
% Q = (sqrt(2 ^ (1 / oct)) / ((2 ^ (1 / oct)) - 1))
% disp('Calculating STFT Spectrogram.')
[x, fs] = loadSignal(7, fftLen);
% [x, fs] = audioread('E:\Download\Sounds of Pine Bark Beetle.wav');
% x = x(1 : 10000);
%% Textbook Morlet CWT
% VoicesPerOctave = 48;
% fb = cwtfilterbank(SignalLength=length(x),Boundary="periodic",WaveletParameters=[2,80],VoicesPerOctave=VoicesPerOctave);
% psif = freqz(fb,FrequencyRange="twosided",IncludeLowpass=true);
% [cfs,~,~,scalcfs] = fb.wt(x);
% imagesc(log10(abs(cfs) + eps))
% xrecAN = ICWT_DualFrame(cfs, psif, scalcfs); % Dual frame form
% xrecSI = icwt(cfs,[],ScalingCoefficients=scalcfs); % Single integral form
% SNR_CWT_FIRLiked_AnalysisSynthesis = 10*log10(sum(abs(x).^2)/sum(abs(x-xrecAN').^2))
% SNR_CWT_SingleIntegral = 10*log10(sum(abs(x).^2)/sum(abs(x-xrecSI').^2))
%% Back to frame-based CWT
% x = [zeros(fftLen + fix(hop / 2), 1); x; zeros(fftLen, 1)];
coeff = coeffCollection{dim1, dim2}.coeff;
f_q = coeffCollection{dim1, dim2}.f_q;
sparCutOff = 1e-10;
gWeights = abs(coeff.correctionMatrix) > sparCutOff;
sparsity = 1 - sum(gWeights(:)) / numel(coeff.correctionMatrix)
disp('Calculating IIR CQT Spectrogram.')
tic
[spec, t_q] = ltv_spectrogram(x, coeff);
% [spec,f_q,t_q,p] = spectrogram([y; zeros(fftLen, 1)],hann(fftLen, 'periodic'),fftLen-hop,fftLen,fs);
toc
% max(20 * log10( abs(spec(:, :, 1)) * 2 ./ (sum(hann(fftLen, 'periodic'))) ), [], 'all')
inverseOvp = 1;
simulateRegularSTFT = 0;
if inverseOvp
    %% Spectral editing
    if simulateRegularSTFT
        editing = spec(:, :, 3) .* 0.25 * (hop / sum(hann(fftLen, 'periodic') .^ 2));
    else
        editing = spec(:, :, 1);
    end
    %     editing(500 : 1200, 120 : 178) = 0;
    %     editedBk = editing;
    if isfield(coeff, 'wndCorrectionWeightingLF') && ~isempty(coeff.correctionMatrix)
        y = ltv_inverse2(editing, coeff);
    elseif isfield(coeff, 'wndCorrectionWeightingLF')
        y = ltv_inverse1(editing, coeff);
    end
    y = y(:);
    y = y(fftLen - hop + 1 : end);
    if length(y) > length(x)
        y = y(1 : length(x));
    else
        x = x(1 : length(y));
    end
    %     [coeff, f_q] = ltv_precomute(1024, 64, fs, oct, order, HFSamplingLimit, 1, 1);
    %     unvealActualMask = ltv_spectrogram(y, coeff);
    %     unvealActualMask=unvealActualMask(:, :, 3);
    %     unvealActualMask = unvealActualMask ./ max(abs(unvealActualMask(:)));
    %     imagesc(t_q, f_q, 20 * log10(abs(unvealActualMask)))
    %     rec = y2 ./ y';
    %     figure(2)
    plot(y);
    %     ylim([0.95, 1.05])
    SNR = 10*log10(sum(abs(x).^2)/sum(abs(x-y).^2)) %calculate the SNR
    % SNR = 20*log10( sqrt(sum(x(:).^2)) / sqrt(sum( (x(:) - y(:)).^2 ))) %calculate the SNR
end
spec2 = spec(:, :, 1);
spec2(halfLen+1:fftLen,:) = conj(spec2(halfLen-1:-1:2,:));
E3 = norm(abs(spec2),'fro').^2;
s_q = spec(:, :, 1) * 2 ./ (fftLen / 2);
% Normalize to 0dB
s_q = s_q ./ max(abs(s_q(:)));
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
imagesc(t_q, f_q, 20 * log10(abs(s_q)))
caxis([-100, 0])
colormap(jet);
set(gca,'YDir','normal');
colorbar
if order == 2
    title('Second order')
else
    title('First order')
end