function spec = ltv_1st_2nd(dftSpec, b, a, c1, c2, prepad, pospad, phaseShifter1, fftLen, halfLen, interpolation)
% Rectangular windowed
x_fft1 = [conj(dftSpec(prepad + 1 : -1 : 2, :)); dftSpec; conj(dftSpec(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
smpShifted = dftSpec .* phaseShifter1;
x_fft2 = [conj(smpShifted(prepad + 1 : -1 : 2, :)); smpShifted; conj(smpShifted(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
%% Gaussian windowing
x_fft = [x_fft1, x_fft2];
tmp = zeros(size(x_fft, 1), 2, 'like', x_fft);
q_fft_frame = ltv(x_fft, tmp, b, a, c1, c2);
% Remove periodic padding
% q_fft_frame(:, 3) = q_fft_frame1;
spec = q_fft_frame(prepad + 1 : end - pospad + 1, :);
% Reference to:
% A high resolution fundamental frequency determination based on phase changes of the Fourier transform
instantaneousFreq = arg(spec(:, 2) .* conj(spec(:, 1))) * fftLen;
instantaneousFreq(instantaneousFreq > halfLen) = fftLen - instantaneousFreq(instantaneousFreq > halfLen);
RS = cplxReassignment(spec(:, 1), instantaneousFreq, interpolation); % Any NaN or Inf will crash
spec(:, 2) = RS;
spec(1, :) = real(spec(1, :));
spec(end, :) = real(spec(end, :));
spec(:, 3) = dftSpec;
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