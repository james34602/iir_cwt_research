wndSize = uint32(linspace(0, 400, 500)) + 1;
st = initCplxMovAvg(wndSize);
sigLen = 39077;
rng(1);
aa = [randn(length(wndSize), 1000), complex(randn(length(wndSize), sigLen), randn(length(wndSize), sigLen)), zeros(length(wndSize), wndSize(end))];
tic
bb = zeros(size(aa), 'like', 1i);
for idx = 1 : sigLen
    bb(:, idx) = procCplxMovAvg(st, aa(:, idx));
end
toc
for idx = 1 : length(wndSize)
    bb(idx, :) = circshift(bb(idx, :), -double(wndSize(idx) / 2 - 1), 2);
end
tic
cc = zeros(size(aa));
for idx = 1 : length(wndSize)
    cc(idx, :) = movmean(aa(idx, :), wndSize(idx));
end
toc
plot(imag(bb(end, :)))
hold on
plot(imag(cc(end, :)))
hold off;
axis tight