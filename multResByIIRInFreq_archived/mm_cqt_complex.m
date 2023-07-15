function [S, f, t] = mm_cqt_complex(x,nfft,hop,fs,Noct)
% transform input data to a column vector
x = x(:);
ny = length(x);
% number of frames (temporal index)
nframes = ceil(ny/hop);
frameindex = 1 + (0 : nframes - 1) * hop; % samples.
% zero padding at the end to complete the last frame.
x = [x; zeros(nfft-mod(ny,hop),1)];
halfLen = nfft/2+1;
% matrix to store the complex CQT-spectrum
S = zeros(nframes, halfLen);
f = (0:1:nfft/2)*fs/nfft;
transformationMatrix = zeros(halfLen, halfLen);
for j = 1:halfLen
    sigma = (f(j) / Noct) / pi; % standard deviation
    if sigma < eps
        sigma = eps;
    end
    g = exp(-( ( (f - f(j) ) .^ 2 ) / (2 * (sigma ^ 2) ) ) ); % Gaussian
    transformationMatrix(j, :) = g ./ sum(g);
end
reconstructionMatrix = inv(transformationMatrix);
wnd = hann(nfft);
% MM CQT transform
for j=1:nframes
    S(j,:) = iir_ltv_q_fft(x(frameindex(j):frameindex(j)+nfft-1) .* wnd, transformationMatrix, reconstructionMatrix);
end
t = (frameindex-1+nfft/2)/fs;
S = S.';
end

function q_fft_frame = iir_ltv_q_fft(x, fb, reconstructionMatrix)
nfft = length(x);
halfLen = nfft/2+1;
x_fft = fft([x(end/2+1:end); x(1:end/2)]);
x_fft = x_fft(1 : halfLen);
q_fft_frame = fb * x_fft;
% rec = reconstructionMatrix * q_fft_frame;
% dif = rec - x_fft;
% error = mean([real(dif), imag(dif)])
end