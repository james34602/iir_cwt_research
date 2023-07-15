wndLen = 512;
smpShift = 128;
wnd = hann(wndLen,'periodic');
kDelta1 = zeros(16384, 1);
pos1 = 997;
len = 2000;
sig = (rand(len, 1) - 0.5) * 2;
kDelta1(pos1 + 1 : pos1 + len) = sig;
spec = mystft(kDelta1, wnd, smpShift, wndLen);
spec(1 : 100, :) = 0;
spec(120 : 180, :) = 0;
spec(220 : end, :) = 0;
rec1 = offlineISTFT(spec, wnd, wndLen - smpShift);

kDelta2 = zeros(16384, 1);
pos2 = 1470;
kDelta2(pos2 + 1 : pos2 + len) = sig;
spec = mystft(kDelta2, wnd, smpShift, wndLen);
spec(1 : 100, :) = 0;
spec(120 : 180, :) = 0;
spec(220 : end, :) = 0;
rec2 = offlineISTFT(spec, wnd, wndLen - smpShift);
rec2 = circshift(rec2, -(pos2 - pos1));
subplot(3,1,1)
plot(rec1);
title('Unshifted')
subplot(3,1,2)
plot(rec2);
title('Shifted')
subplot(3,1,3)
plot(rec1 - rec2);
title('Difference')
mae = sum((rec1 - rec2) .* (rec1 - rec2))
function [S, frontPad, endPad] = mystft(x, window, hop, NFFT)
%% default values
if nargin<2 || isempty(window)
    window = hann(128,'periodic');
end
windowLen = length(window);

if nargin<3 || isempty(hop)
    hop = max(1, round(windowLen/4));
end

if nargin<4 || isempty(NFFT)
    NFFT = 2^nextpow2(windowLen);
end

%% Padd with zeros from the left (for causality)
x = [zeros(windowLen / 2, 1); x;];
%% partition x to data segments
n = size(x,1);
I = 1:hop:n;
frontPad = windowLen / 2;
endPad = (I(end)+windowLen-1) - n;
m = (0:windowLen-1)';
I = I+m;
x(end+1:I(end), :)=0;
x = x(I);

%% apply window function
x = x.*window(:);

%% FFT each segment
S = fft(x, NFFT, 1);
S = S(1 : NFFT / 2 + 1, :);
end
function [output, inverseTransformStability]=offlineISTFT(targetSpectrogram,window,noverlap)
N=(size(targetSpectrogram,1)-1)*2;
hop=N-noverlap;

frequencyframe=zeros(N,1,class(targetSpectrogram));

nframes=size(targetSpectrogram,2);
output=zeros(nframes * hop + (N - hop), 1);

%% Least-squares correction
L2 = ceil(N/hop)*hop;
window2 = window;
window2(end+1:L2) = 0; % zero pad window to be an integer multiple of hop
J = mod( (1:hop:L2)' + (0:hop-1) -1, L2 ) + 1;
denom = 1.0 ./ (sum(window2(J).^2, 1)');
inverseTransformStability = max(denom) / min(denom);
t=1;
for i=1:nframes
	frequencyframe(1:N/2+1) = targetSpectrogram(1:N/2+1,i);
    frequencyframe(N/2+2:N) = flipud(real(frequencyframe(2:N/2)) -1j * imag(frequencyframe(2:N/2)));
    
    timeframe = real(ifft(frequencyframe));
    framed = window .* timeframe;
    output(t:t+N-1) = output(t:t+N-1) + framed;
	t = t + hop;
end
denom = repmat(denom, ceil((nframes * hop + (N - hop))/hop), 1);
denom((nframes * hop + (N - hop))+1:end) = [];
output = output .* denom;
end