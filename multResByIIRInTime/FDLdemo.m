wndSize = linspace(0, 50, 300);
st = fractionalDL(wndSize);
sigLen = 7700;
rng(1);
aa = [randn(length(wndSize), 10), complex(randn(length(wndSize), sigLen), randn(length(wndSize), sigLen)), zeros(length(wndSize), wndSize(end))];
tic
bb = zeros(size(aa), 'like', 1i);
for idx = 1 : length(aa)
    bb(:, idx) = st.process(aa(:, idx));
end
tic
cc = zeros(size(aa), 'like', 1i);
for idx = 1 : length(wndSize)
    cc(idx, :) = circshift(aa(idx, :), fix(wndSize(idx)), 2);
end
toc
idx = 8;
disp("Delay = " + string(wndSize(idx)))
plot(imag(bb(idx, :)))
hold on
plot(imag(cc(idx, :)))
hold off;
axis tight