lennon = 128;
order = 4;
coeffB = zeros(lennon, 1);
coeffA = zeros(order + 1, lennon);
for idx = 1 : lennon
    [coeffB(idx, :), coeffA(:, idx)] = bessel(order, idx);
end
st = initNPole1Zero(uint32(lennon), uint32(order));
sigLen = 40700;
rng(1);
aa = [ones(lennon, 1), complex(randn(lennon, sigLen), rand(lennon, sigLen)), complex(ones(lennon, 1), -ones(lennon, 1) * 0.4), zeros(lennon, 800)];
tic
bb = zeros(size(aa), 'like', 1i);
for idx = 1 : length(aa)
    bb(:, idx) = procNPole1Zero(st, coeffB, coeffA, aa(:, idx));
end
toc
fvtool(bb(1, :))
fvtool(bb(end/4, :))
fvtool(bb(end, :))