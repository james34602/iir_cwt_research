function spec = ltv_1st_2nd(dftSpec, b, a, c1, c2, prepad, pospad, phaseShifter1, phaseShifter2, phaseShifter3, corrF, fftLen, halfLen, interpolation)
%% Hann in frequency domain with input being shifted by fftshift
specHann = 2 * dftSpec + [conj(dftSpec(2)); dftSpec(1 : end - 1)] + [dftSpec(2 : end); conj(dftSpec(end - 1))];
% specHann = 2 * specHann + [conj(specHann(2)); specHann(1 : end - 1)] + [specHann(2 : end); conj(specHann(end - 1))];
%% Hann in frequency domain with 1 sample delayed input being shifted by fftshift
shiftedspecHann = phaseShifter1 .* (2 * dftSpec + phaseShifter2 .* [conj(dftSpec(2)); dftSpec(1 : end - 1)] + phaseShifter3 .* [dftSpec(2 : end); conj(dftSpec(end - 1))]);
% Hann windowed
x_fft1 = [conj(specHann(prepad + 1 : -1 : 2, :)); specHann; conj(specHann(halfLen - 1 : -1 : (halfLen - pospad + 1), :))] / 4;
x_fft2 = [conj(shiftedspecHann(prepad + 1 : -1 : 2, :)); shiftedspecHann; conj(shiftedspecHann(halfLen - 1 : -1 : (halfLen - pospad + 1), :))] / 4;
% Rectangular windowed
% x_fft1 = [conj(dftSpec(prepad + 1 : -1 : 2, :)); dftSpec; conj(dftSpec(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
% smpShifted = dftSpec .* phaseShifter1;
% x_fft2 = [conj(smpShifted(prepad + 1 : -1 : 2, :)); smpShifted; conj(smpShifted(halfLen - 1 : -1 : (halfLen - pospad + 1), :))];
%% Gaussian windowing
x_fft = [x_fft1, x_fft2];
tmp = zeros(size(x_fft, 1), 2, 'like', x_fft);
q_fft_frame = ltv(x_fft, tmp, b, a, c1, c2);
% Remove periodic padding
spec = q_fft_frame(prepad + 1 : end - pospad + 1, :);
% Reference to:
% A high resolution fundamental frequency determination based on phase changes of the Fourier transform
instantaneousFreq = arg(spec(:, 2) .* conj(spec(:, 1))) * fftLen;
instantaneousFreq(instantaneousFreq > halfLen) = fftLen - instantaneousFreq(instantaneousFreq > halfLen);
RS = cplxReassignment(spec(:, 1), instantaneousFreq, interpolation); % Any NaN or Inf will crash
spec(:, 2) = RS;
if ~isempty(corrF)
    spec = spec .* corrF;
end
end
function y = arg(values)
% Argument (angle) of complex numbers wrapped and scaled to [0.0, 1.0].
y = mod(angle(values) / (2 * pi), 1.0);
end