function y = overlapAdd(tmp, hop)
nframes = size(tmp, 2);
fftLen = size(tmp, 1);
xlen = fftLen + (nframes-1)*hop;
y = zeros(xlen, 1);
for l = 1 : nframes
    y(1+(l-1)*hop : fftLen+(l-1)*hop) = y(1+(l-1)*hop : fftLen+(l-1)*hop) + tmp(:, l);
end
end