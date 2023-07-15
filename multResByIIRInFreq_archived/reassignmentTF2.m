load ../signal1
x = s;
% x(:)=0;
% x(:)=1;
reassign_time = 0;
reassign_frequency = 1;
tfrMag = 0;
x = [x; zeros(32, 1)];
sigLen = length(x);
fftLen = 64;
halfLen = fftLen / 2 + 1;
hopSize = 16;
x_frames = buffer(x, fftLen, fftLen - hopSize);
start_times = linspace(1, size(x_frames, 2), ceil(sigLen / hopSize));
w = hann(fftLen);
% x_frames(:)=1;
unwnd_complex = fft(x_frames);
% spectra shifted in frequency by +1 spectral linr
X_prev_freq_p1 = circshift(unwnd_complex, 1);
% spectra shifted in frequency by -1 spectral line
X_prev_freq_n1 = circshift(unwnd_complex, -1);
STFT_th = X_prev_freq_p1 - X_prev_freq_n1;
% recTmd = imag(ifft(STFT_th));
% plot(recTmd(:, 26))
X_complex = fft(x_frames .* w);

fout=((0:(fftLen-1))/fftLen*2*pi)';
fmax = fout(end);
fcorr2 = -imag(((1j * -STFT_th)*0.5*0.049768449558477) ./ X_complex);
fcorr2(~isfinite(fcorr2)) = 0;
fcorr2 = fout + fcorr2;
rearrangement2 = 1 + mod(round(fcorr2 * (fftLen - 1) / fmax), fftLen);
dwnd = diff(w);
%%
%% Derivative of main window
first = [dwnd(1) ; dwnd];
last = [dwnd; dwnd(end)];
dwindow = ((first + last) / 2);
dX_complex = fft(x_frames .* dwindow);
%% SST
fcorr = -imag(dX_complex ./ X_complex);
% fcorr = -imag((dX_complex .* conj(X_complex)) ./ (abs(X_complex) .^ 2));
fcorr(~isfinite(fcorr)) = 0;
fcorr = fout + fcorr;
% Sxx = X1 .* ez;
% tfr(:, nFrame) = Sxx(1 : fftLen / 2 + 1);
rearrangement = 1 + mod(round(fcorr * (fftLen - 1) / fmax), fftLen);
%%
% spectra shifted in frequency
X_prev_freq = circshift(X_complex, 1);
% spectra of signal shifted in time
% This fakes looking at the previous frame shifted by one sample
cirX_frames = circshift(x_frames, 1); % Fake delay
X_prev_time = fft(cirX_frames .* w);
% X_complex = X_complex(1 : halfLen, :);
% X_prev_time = X_prev_time(1 : halfLen, :);
% X_prev_freq = X_prev_freq(1 : halfLen, :);
X_mag = abs(X_complex);
% cross-spectra - ie. spectra of cross-correlation between the
% respective time-domain signals
X_cross_time = cross_spectrum(X_complex, X_prev_time);
X_cross_freq = cross_spectrum(X_complex, X_prev_freq);
% instantaneous frequency estimates
% normalized frequencies in range [0.0, 1.0] - from DC to sample rate
X_inst_freqs = arg(X_cross_time);
X_inst_freqs2 = arg(X_prev_time) - arg(X_complex);
X_inst_freqs2(X_inst_freqs2 < 0) = 1 + X_inst_freqs2(X_inst_freqs2 < 0);
X_inst_freqs2 = 1 + -X_inst_freqs2;
% X_inst_freqs = X_inst_freqs2;
% plot(X_inst_freqs(:, 11))
% hold on
% plot(X_inst_freqs2(:, 11));
% hold off
% axis tight
% instantaneous group delay estimates
% relative coordinates within the frame with range [-0.5, 0.5] where
% 0.0 is the frame center
X_group_delays = 0.5 - arg(X_cross_freq);
X_group_delays2 = arg(X_prev_freq) - arg(X_complex);
% X_group_delays2(X_group_delays2 > 0.5) = 1 - X_group_delays2(X_group_delays2 > 0.5);
X_group_delays2(X_group_delays2 < 0) = 1 + X_group_delays2(X_group_delays2 < 0);
X_group_delays2 = -0.5 + X_group_delays2;
% X_group_delays2 = X_group_delays2 / fftLen;
% plot(X_group_delays(:, 265))
% hold on
% plot(X_group_delays2(:, 265));
% hold off
% axis tight
% Reassigned spectrogram requantized both in frequency and time.
% Note it is quantized into non-overlapping output time frames which may be
% of a different size than input time frames.
% group delays are in range [-0.5, 0.5] - relative coordinates within the
% frame where 0.0 is the frame center
X_time = repmat(start_times, [fftLen, 1]);
if reassign_time
    X_time = X_time + X_group_delays * (fftLen / hopSize);
end
if reassign_frequency
    X_y = X_inst_freqs;
else
    X_y = repmat(linspace(0, 1, fftLen), [size(X_time, 2), 1]).';
end
inds = 0 : (fftLen - 1);
m = floor(fftLen / 2);
ez = exp(-1i*2*pi*m*inds/fftLen)';
X_complex = X_complex .* ez;
X_time = X_time';
X_y = X_y';
X_mag = X_mag';
X_complex = X_complex.';
sz = fliplr(size(X_time));
if tfrMag
    Snew = X_mag;
else
    Snew = X_complex;
end
RS = zeros(sz);
range1 = round(X_time);
range1(range1 < 1) = 1;
range1(range1 > size(X_time, 1)) = size(X_time, 1);
range2 = 1 + mod(round(X_y * size(X_time, 2)), fftLen);
range2 = rearrangement2';
for k = 1 : size(X_time, 1)
    mag = Snew(k, :);
    for j = 1 : fftLen
        RS(range2(k, j), range1(k, j)) = RS(range2(k, j), range1(k, j)) + mag(j);
    end
end
if tfrMag
    imagesc(log10(abs(RS) + 0.000004))
else
    imagesc(log10(abs(RS) + 0.001))
end
colormap(jet);
set(gca,'YDir','normal');
function y = cross_spectrum(spectrumA, spectrumB)
%     Returns a cross-spectrum, ie. spectrum of cross-correlation of two signals.
%     This result does not depend on the order of the arguments.
%     Since we already have the spectra of signals A and B and and want the
%     spectrum of their cross-correlation, we can replace convolution in time
%     domain with multiplication in frequency domain.
y = spectrumA .* conj(spectrumB);
end
function y = arg(values)
%     Argument (angle) of complex numbers wrapped and scaled to [0.0, 1.0].
%
%     input: an array of complex numbers
%     output: an array of real numbers of the same shape
%
%     np.angle() returns values in range [-np.pi, np.pi].
phNorm = angle(values) / (2 * pi);
% y = mod(angle(values) / (2 * pi), 1.0);
y = phNorm;
y(y < 0) = (y(y < 0) + 1);
end