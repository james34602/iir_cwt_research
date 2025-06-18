function y = overlapAdd3D(tmp, hop)
tmp=permute(tmp, [1, 3, 2]);
nframes = size(tmp, 3);
nchannels = size(tmp, 2);
fftLen = size(tmp, 1);
xlen = fftLen + (nframes-1)*hop;
y = zeros(xlen, nchannels);
for l = 1 : nframes
    y(1+(l-1)*hop : fftLen+(l-1)*hop, :) = y(1+(l-1)*hop : fftLen+(l-1)*hop, :) + tmp(:, :, l);
end
end