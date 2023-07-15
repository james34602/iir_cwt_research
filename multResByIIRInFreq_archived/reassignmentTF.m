load ../signal1
x = s;
reassign_time = 1;
reassign_frequency = 1;
fakeDelay = 0;
tfrMag = 0;
x = [x; zeros(32, 1)];
sigLen = length(x);
frameSize = 64;
halfLen = frameSize / 2 + 1;
hopSize = 1;
x_frames = buffer(x, frameSize, frameSize - hopSize);
start_times = linspace(1, size(x_frames, 2), ceil(sigLen / hopSize));
w = hann(frameSize);
X_complex = fft(x_frames .* w);
% linear magnitude spectrogram
X_mag = abs(X_complex) / frameSize;
% spectra shifted in frequency
X_prev_freq = circshift(X_complex, 1);
% spectra of signal shifted in time
% This fakes looking at the previous frame shifted by one sample
if fakeDelay
    cirX_frames = circshift(x_frames, 1); % Fake delay
    X_prev_time = fft(cirX_frames .* w);
else
    delayed = circshift(x, 1);
    delayed(1) = 0;
    delayedX_frames = buffer(delayed, frameSize, frameSize - hopSize);
    X_prev_time = fft(delayedX_frames .* w);
end
% cross-spectra - ie. spectra of cross-correlation between the
% respective time-domain signals
X_cross_time = cross_spectrum(X_complex, X_prev_time);
X_cross_freq = cross_spectrum(X_complex, X_prev_freq);
% instantaneous frequency estimates
% normalized frequencies in range [0.0, 1.0] - from DC to sample rate
X_inst_freqs = arg(X_cross_time);
% instantaneous group delay estimates
% relative coordinates within the frame with range [-0.5, 0.5] where
% 0.0 is the frame center
X_group_delays = 0.5 - arg(X_cross_freq);

X_mag = X_mag(1 : halfLen, :);
% X_mag(1, :) = X_mag(1, :) * 0.5;
X_complex = X_complex(1 : halfLen, :);
X_inst_freqs = X_inst_freqs(1 : halfLen, :);
X_group_delays = X_group_delays(1 : halfLen, :);
% Reassigned spectrogram requantized both in frequency and time.
% Note it is quantized into non-overlapping output time frames which may be
% of a different size than input time frames.
% group delays are in range [-0.5, 0.5] - relative coordinates within the
% frame where 0.0 is the frame center
X_time = repmat(start_times, [halfLen, 1]);
if reassign_time
    X_time = X_time + X_group_delays * (frameSize / hopSize);
end
if reassign_frequency
    X_y = X_inst_freqs;
else
    X_y = repmat(linspace(0, 0.5, halfLen), [size(X_time, 2), 1]).';
end
inds = 0 : (frameSize - 1);
m = floor(frameSize / 2);
ez = exp(-1i*2*pi*m*inds/frameSize)';
X_complex = X_complex .* ez(1 : halfLen);
X_time = X_time';
X_y = X_y';
X_mag = X_mag';
X_complex = X_complex.';
That = X_time(:) / max(X_time(:)) * size(X_time, 1);
That(That < 1) = 1;
That(That > size(X_time, 1)) = size(X_time, 1);
Fhat = X_y(:) * size(X_time, 2) * 2 + 1;
Fhat(Fhat < 1) = 1;
Fhat(Fhat > size(X_time, 2)) = size(X_time, 2);
sz = fliplr(size(X_time));
if tfrMag
    Snew = X_mag(:);
else
    Snew = X_complex(:);
end
% RS = zeros(sz);
% range1 = round(X_time);
% range1(range1 < 1) = 1;
% range1(range1 > size(X_time, 1)) = size(X_time, 1);
% range2 = round(X_y * size(X_time, 2) * 2) + 1;
% range2(range2 < 1) = 1;
% range2(range2 > size(X_time, 2)) = size(X_time, 2);
% for k = 1 : size(X_time, 1)
%     mag = X_mag(k, :);
%     for j = 1 : halfLen
%         RS(range2(k, j), range1(k, j)) = RS(range2(k, j), range1(k, j)) + mag(j);
%     end
% end
RS = accumarray([round(Fhat), round(That)],Snew,sz);
if tfrMag
    imagesc(log10(abs(RS) + 0.000004))
else
    imagesc(log10(abs(RS) + 0.01))
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
y = mod(angle(values) / (2 * pi), 1.0);
end