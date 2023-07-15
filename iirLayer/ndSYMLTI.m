rng(1)
clear
m = 4;
nSigs = 3;
input = randn(m, nSigs);
dp = randn(nSigs, 1);
nSections = 2;
coefNum = 1;
parallelFilter = 1;
reverseFilter = 1;
statesNoGrad = 0;
target = randn(m, 1);
%% My implementation
b0 = rand(coefNum, nSections);
d1 = rand(coefNum, nSections);
d2 = rand(coefNum, nSections);
c1 = rand(coefNum, nSections);
c2 = rand(coefNum, nSections);
initialValueZ1 = rand(nSigs, nSections);
initialValueZ2 = rand(nSigs, nSections);
maspoles = [input(:); b0(:); d1(:); d2(:); c1(:); c2(:)];
if ~statesNoGrad
    maspoles = [maspoles; initialValueZ1(:); initialValueZ2(:)];
end
poles = sym('poles', [length(maspoles), 1]);
flat = poles(5 * (nSections - 1) + 5 + 1 : 5 * (nSections - 1) + 5 + m * nSigs);
input = reshape(poles(1 : m * nSigs), [m, nSigs]);
b0 = reshape(poles(m * nSigs + coefNum * nSections * 0 + 1 : m * nSigs + coefNum * nSections * 1), [coefNum, nSections]);
d1 = reshape(poles(m * nSigs + coefNum * nSections * 1 + 1 : m * nSigs + coefNum * nSections * 2), [coefNum, nSections]);
d2 = reshape(poles(m * nSigs + coefNum * nSections * 2 + 1 : m * nSigs + coefNum * nSections * 3), [coefNum, nSections]);
c1 = reshape(poles(m * nSigs + coefNum * nSections * 3 + 1 : m * nSigs + coefNum * nSections * 4), [coefNum, nSections]);
c2 = reshape(poles(m * nSigs + coefNum * nSections * 4 + 1 : m * nSigs + coefNum * nSections * 5), [coefNum, nSections]);
states = cell(nSections, 1);
if ~statesNoGrad
    initialValueZ1 = reshape(poles(m * nSigs + coefNum * nSections * 5 + 1 : m * nSigs + coefNum * nSections * 5 + nSigs * nSections), [nSigs, nSections]);
    initialValueZ2 = reshape(poles(m * nSigs + coefNum * nSections * 5 + nSigs * nSections + 1 : m * nSigs + coefNum * nSections * 5 + nSigs * nSections * 2), [nSigs, nSections]);
    for idx = 1 : nSections
        states{idx}.z1_A = initialValueZ1(:, idx).';
        states{idx}.z2_A = initialValueZ2(:, idx).';
    end
else
    initialValueZ1 = 0;
    initialValueZ2 = 0;
    for idx = 1 : nSections
        states{idx}.z1_A = zeros(1, nSigs);
        states{idx}.z2_A = zeros(1, nSigs);
    end
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
sub = tmp2 * dp - target;
mse1 = sum(sub.^2);
gradf = jacobian(mse1, poles);
fh = matlabFunction(mse1,gradf,'vars',{poles});
[mse1, gradf] = fh(maspoles);
plot(gradf);axis tight