rng(1)
clear
m = 4;
nSigs = 3;
ltv = 0;
if ltv == 1
    coefNum = m;
else
    coefNum = 1;
end
parallelFilter = 0;
reverseFilter = 0;
[tarB, tarA] = butter(3,[0.2 0.6],'stop');
sosIdeal = tf2sos(tarB, tarA);
nSections = size(sosIdeal, 1);
input = randn(m, nSigs);
target = filter(tarB, tarA, input);
%% My implementation
[tarB, tarA] = cheby1(3,5, [0.2 0.6],'stop');
sos = tf2sos(tarB, tarA);
b0 = sos(:, 1).';
b1 = sos(:, 2).';
b2 = sos(:, 3).';
a1 = sos(:, 5).';
a2 = sos(:, 6).';
c1 = a1 + 2;
c2 = (1 + a1 + a2) ./ c1;
d1 = (2 * b0 + b1) ./ c1;
d2 = (b0 + b1 + b2) ./ (c1 .* c2);
maspoles = [b0(:); d1(:); d2(:); c1(:); c2(:)];
%%
poles = sym('poles', [length(maspoles), 1]);
b0 = reshape(poles(coefNum * nSections * 0 + 1 : coefNum * nSections * 1), [coefNum, nSections]);
d1 = reshape(poles(coefNum * nSections * 1 + 1 : coefNum * nSections * 2), [coefNum, nSections]);
d2 = reshape(poles(coefNum * nSections * 2 + 1 : coefNum * nSections * 3), [coefNum, nSections]);
c1 = reshape(poles(coefNum * nSections * 3 + 1 : coefNum * nSections * 4), [coefNum, nSections]);
c2 = reshape(poles(coefNum * nSections * 4 + 1 : coefNum * nSections * 5), [coefNum, nSections]);
states = cell(nSections, 1);
initialValueZ1 = 0;
initialValueZ2 = 0;
for idx = 1 : nSections
    states{idx}.z1_A = zeros(1, nSigs);
    states{idx}.z2_A = zeros(1, nSigs);
end
% Forward filtering
if parallelFilter == 1
    if reverseFilter
        for a = m : -1 : 2
            in = input(a, :);
            for idx = 1 : nSections
                y1 = in - states{idx}.z1_A - states{idx}.z2_A;
                out(idx, :) = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
                states{idx}.z2_A = states{idx}.z2_A + c2(idx) .* states{idx}.z1_A;
                states{idx}.z1_A = states{idx}.z1_A + c1(idx) .* y1;
            end
            tmp2(a, :) = sum(out, 1);
        end
        in = input(1, :);
        for idx = 1 : nSections
            y1 = in - states{idx}.z1_A - states{idx}.z2_A;
            out(idx, :) = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
        end
        tmp2(1, :) = sum(out, 1);
    else
        for a = 1 : m - 1
            in = input(a, :);
            for idx = 1 : nSections
                y1 = in - states{idx}.z1_A - states{idx}.z2_A;
                out(idx, :) = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
                states{idx}.z2_A = states{idx}.z2_A + c2(idx) .* states{idx}.z1_A;
                states{idx}.z1_A = states{idx}.z1_A + c1(idx) .* y1;
            end
            tmp2(a, :) = sum(out, 1);
        end
        in = input(m, :);
        for idx = 1 : nSections
            y1 = in - states{idx}.z1_A - states{idx}.z2_A;
            out(idx, :) = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
        end
        tmp2(m, :) = sum(out, 1);
    end
else
    if reverseFilter
        for a = m : -1 : 2
            out = input(a, :);
            for idx = 1 : nSections
                y1 = out - states{idx}.z1_A - states{idx}.z2_A;
                out = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
                states{idx}.z2_A = states{idx}.z2_A + c2(idx) .* states{idx}.z1_A;
                states{idx}.z1_A = states{idx}.z1_A + c1(idx) .* y1;
            end
            tmp2(a, :) = out;
        end
        out = input(1, :);
        for idx = 1 : nSections
            y1 = out - states{idx}.z1_A - states{idx}.z2_A;
            out = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
        end
        tmp2(1, :) = out;
    else
        for a = 1 : m - 1
            out = input(a, :);
            for idx = 1 : nSections
                y1 = out - states{idx}.z1_A - states{idx}.z2_A;
                out = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
                states{idx}.z2_A = states{idx}.z2_A + c2(idx) .* states{idx}.z1_A;
                states{idx}.z1_A = states{idx}.z1_A + c1(idx) .* y1;
            end
            tmp2(a, :) = out;
        end
        out = input(m, :);
        for idx = 1 : nSections
            y1 = out - states{idx}.z1_A - states{idx}.z2_A;
            out = b0(idx) .* y1 + d1(idx) .* states{idx}.z1_A + d2(idx) .* states{idx}.z2_A;
        end
        tmp2(m, :) = out;
    end
end
sub = tmp2 - target;
mse1 = mean(sub(:).^2);
gradf = jacobian(mse1, poles);
fh = matlabFunction(mse1,gradf,'vars',{poles});
[mse1, gradf] = fh(maspoles);