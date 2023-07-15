fftLen = 4096;
hop = 256;
% [y, fs] = audioread('../sinesweep.wav');
% 
% disp('Calculating STFT Spectrogram.')
% [s,f,t,p] = spectrogram([y; zeros(fftLen, 1)],hann(fftLen),fftLen-hop,fftLen,fs); 

Noct = 28;
disp('Calculating IIR CQT Spectrogram.')
tic
[s_mmq, f_mmq, t_mmq] = mm_cqt_complex(y,fftLen,hop,fs,Noct);
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

function [S, f, t] = mm_cqt_complex(x,fftLen,hop,fs,Noct)
% transform input data to a column vector
x = x(:);
xBk = x;
x = [zeros(fftLen - hop, 1); x; zeros(fftLen / 2, 1)];
ny = length(x);
% number of frames (temporal index)
nframes = ceil(ny/hop);
frameindex = 1 + (0 : nframes - 1) * hop; % samples.
% zero padding at the end to complete the last frame.
x = [x; zeros(fftLen-mod(ny,hop),1)];
halfLen = fftLen/2+1;
% matrix to store the complex CQT-spectrum
S = zeros(nframes, halfLen);
f = (0:1:fftLen/2)*fs/fftLen;
% number of points of pre and post padding used to set initial conditions
prepad = 100;
pospad = 100;
% number of points of pre and post padding used to set initial conditions
thetas1 = (0:(fftLen/2+pospad)-1);
thetas1(1) = eps;
thetas1 = [fliplr(thetas1(2 : prepad + 1)), thetas1];
thetas1 = thetas1(1:fftLen/2+prepad+1);
thetas1 = [thetas1, thetas1(length(thetas1)-1:-1:1)];
thetas1 = thetas1(1 : halfLen + prepad + pospad - 1);
thetas1 = thetas1 ./ thetas1(halfLen+prepad);
thetas1 = thetas1 * fs / 2;
%%
paddedSpecLen = length(thetas1);
sigma = (thetas1 ./ Noct) / pi; % standard deviation
sigma(sigma > 100) = 100;
[bSm, aSm] = butter(1, 0.1);
sigma = filtfilt(bSm, aSm, sigma);
% plot(sigma)
transformationMatrix = zeros(paddedSpecLen, paddedSpecLen);
transformationMatrix2 = zeros(paddedSpecLen, paddedSpecLen);
ovp = fftLen / hop;
for j = 1:paddedSpecLen
%     if j < prepad + fftLen/4
%         sigma(j) = 100;
%     else
%         sigma(j) = 150;
%     end
    g = exp(-( ( (thetas1 - thetas1(fftLen / 4)  ) .^ 2 ) / (2 * (sigma(j) ^ 2) ) ) ); % Gaussian
    firstValid = find(g > 0.0, 1, 'first');
    lastValid = find(g > 0.0, 1, 'last');
    validLen = lastValid - firstValid;
%     plot(g(firstValid : lastValid))
%     g = g ./ sum(g);
    %% Generate correction window
    gaussFilt = fft(g, fftLen);
    gTime = abs(fftshift(gaussFilt));
    g = circshift(g, j - fftLen / 4);
    if (j - fftLen / 4) < fftLen / 4
        g(j + validLen + 1 : end) = 0;
    end
    if (j - fftLen / 4) > fftLen / 4
        g(1 : end - j + validLen) = 0;
    end
    %%
    correctionWndHF = overlapAdd(repmat(gTime, ovp * 2, 1 )', hop);
    correctionWndHF = correctionWndHF(fftLen - hop + 1 : fftLen * 2 - hop);
    correctionWndHF = 1 ./ correctionWndHF;
    halfCorrWnd = correctionWndHF(1 : halfLen);
    halfCorrWnd(halfLen+1:fftLen) = conj(halfCorrWnd(halfLen-1:-1:2));
%     halfCorrWnd = correctionWndHF;
    invHalfCorr = ifft(halfCorrWnd);
    invHalfCorr = ifftshift(invHalfCorr);
    invHalfCorr(abs(invHalfCorr) < 1e-14) = 0;
    g2 = circshift(invHalfCorr, j - fftLen / 2 - 1);
    g2 = g2(1 : prepad + pospad + halfLen - 1);
    if (j - fftLen / 4) < fftLen / 4
        g2(j + validLen + 1 : end) = 0;
    end
    if (j - fftLen / 4) > fftLen / 4
        g2(1 : end - j + validLen) = 0;
    end
    transformationMatrix(:, j) = g;
    transformationMatrix2(:, j) = g2;
%     if j > 2000
%         plot(g2);
%         axis tight;
%         disp('')
%     end
end
% err = Inf;
% transformationMatrix3 = zeros(paddedSpecLen, paddedSpecLen);
% while err > 1e-4
%     transformationMatrix4 = transformationMatrix * transformationMatrix2;
%     err = 0;
%     for j = 1:paddedSpecLen
%         filter800 = transformationMatrix4(:, j);
%         gaussFilt = fft(filter800, fftLen);
%         gTime = abs(fftshift(gaussFilt))';
%         recon = overlapAdd(repmat(gTime, ovp * 2, 1 )', hop);
%         recon = recon(fftLen - hop + 1 : fftLen * 2 - hop);
%         err = err + (max(recon) - min(recon));
%         recon = 1 ./ recon;
%         halfCorrWnd = recon(1 : halfLen);
%         halfCorrWnd(halfLen+1:fftLen) = conj(halfCorrWnd(halfLen-1:-1:2));
% %         plot(recon)
%         idealCorred = ifft(halfCorrWnd);
%         idealCorred = ifftshift(idealCorred);
%         g3 = circshift(idealCorred, j - fftLen / 2 - 1);
%         g3 = g3(1 : prepad + pospad + halfLen - 1);
%         if (j - fftLen / 4) < fftLen / 4
%             g3(j + validLen + 1 : end) = 0;
%         end
%         if (j - fftLen / 4) > fftLen / 4
%             g3(1 : end - j + validLen) = 0;
%         end
%         transformationMatrix3(:, j) = g3;
%     end
%     err = err / paddedSpecLen;
%     disp(err)
%     transformationMatrix2 = transformationMatrix2 * transformationMatrix3;
%     disp('')
% end
% x_fft = fft(ones(fftLen, 1));
% x_fft = x_fft(1 : halfLen);
% spec = iir_ltv_q_fft(x_fft, transformationMatrix, prepad, pospad, transformationMatrix2);
% q_fft_frameinv = spec.';
% q_fft_frameinv(1) = real(q_fft_frameinv(1));
% q_fft_frameinv(end) = real(q_fft_frameinv(end));
% q_fft_frameinv(halfLen+1:fftLen) = conj(q_fft_frameinv(halfLen-1:-1:2));
% recoveredFrame = ifftshift(ifft(q_fft_frameinv));
% plot(recoveredFrame);hold on;plot(gTime .* correctionWndHF');hold off;
% MM CQT transform
transformationMatrix = transformationMatrix * transformationMatrix2;
accFrame2 = zeros(fftLen, 1);
y = zeros(nframes * hop, 1);
for j=1:nframes
    curFrame = x(frameindex(j):frameindex(j)+fftLen-1);
    x_fft = fft(fftshift(curFrame));
    x_fft = x_fft(1 : halfLen);
    S(j,:) = iir_ltv_q_fft(x_fft, transformationMatrix, prepad, pospad, []);
    q_fft_frameinv = S(j,:).';
    q_fft_frameinv(1) = real(q_fft_frameinv(1));
    q_fft_frameinv(end) = real(q_fft_frameinv(end));
    q_fft_frameinv(halfLen+1:fftLen) = conj(q_fft_frameinv(halfLen-1:-1:2));
    recoveredFrame = ifftshift(ifft(q_fft_frameinv));
%     plot(recoveredFrame);axis tight
    % Overlap add
    accFrame2 = accFrame2 + recoveredFrame;
    myOut = accFrame2(1 : hop);
    accFrame2 = [accFrame2(hop + 1 : end); zeros(hop, 1)];
    y(frameindex(j):frameindex(j)+hop-1) = myOut;
end
y = y(fftLen - hop + 1 : end - fftLen / 2);
t = (frameindex-1+fftLen/2)/fs;
S = S.';
end

function [spec, q_fft_frame] = iir_ltv_q_fft(x_fft, fb, prepad, pospad, fb2)
halfLen = length(x_fft);
x_fft = [conj(x_fft(prepad + 1 : -1 : 2, :)); x_fft; conj(x_fft(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
if ~isempty(fb)
    q_fft_frame = fb * x_fft;
else
    q_fft_frame = x_fft;
end
if ~isempty(fb2)
%     q_fft_frame = q_fft_frame(prepad + 1 : end - pospad + 1, :);
%     q_fft_frame = [conj(q_fft_frame(prepad + 1 : -1 : 2, :)); q_fft_frame; conj(q_fft_frame(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
    q_fft_frame = fb2 * q_fft_frame;
end
spec = q_fft_frame(prepad + 1 : end - pospad + 1, :);

end
function [mse, grad] = fitCosine(n, idx, w, target)
len = length(idx);
recipocal = 1 ./ len;
d = 0.5 / (n - 1);
x = idx * d;
err = target - (0.5 - 0.5 * cos(2 * pi * x));
mse = mean(err .* err .* w);
if nargout > 1
    grad = 2*sum((recipocal .* idx .* 0.5 .* pi .* w .* sin((idx .* pi) ./ (n - 1)) .* (target - 0.5 + 0.5 .* cos((idx .* pi) ./ (n - 1)))) ./ (n - 1).^2);
end
end