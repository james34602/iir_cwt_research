rng(1)
fs = 10e3;
wndLen = 256;
x = [zeros(1, wndLen), zeros(1, wndLen * 4)];
x(:) = 0;
x(wndLen + 1) = 1;
x(:)=randn(size(x));
wnd = hann(wndLen, 'periodic');
hop1 = wndLen - wndLen / 64;
hop2 = wndLen - wndLen / 2;
[spec1,f,t1] = stft(x,fs,Window=wnd,OverlapLength=hop1,FFTLength=wndLen);
[spec2,f,t2] = stft(x,fs,Window=wnd,OverlapLength=hop2,FFTLength=wndLen);
spec1Mag = abs(spec1);
spec2Mag = abs(spec2);
gm1 = spec1Mag' * spec1Mag;
gm2 = spec2Mag' * spec2Mag;
function y = overlapAdd(tmp, hop)
nframes = size(tmp, 2);
fftLen = size(tmp, 1);
xlen = fftLen + (nframes-1)*hop;
y = zeros(xlen, 1);
for l = 1 : nframes
    y(1+(l-1)*hop : fftLen+(l-1)*hop) = y(1+(l-1)*hop : fftLen+(l-1)*hop) + tmp(:, l);
end
end