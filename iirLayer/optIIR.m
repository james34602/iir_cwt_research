% rng(1)
clear
m = 1000;
ltv = 0;
if ltv == 1
    coefNum = m;
else
    coefNum = 1;
end
parallelFilter = 0;
reverseFilter = 0;
[tarB, tarA] = butter(3,[0.02 0.1],'stop');
sosIdeal = tf2sos(tarB, tarA);
nSections = size(sosIdeal, 1);
% input = randn(m, nSigs);
fs=48000;
halfLen = 513;fftLen=1024;
cpLen=m;D = cpLen / fs;t = (0:cpLen-1) / fs;f1 = 1 * fs / fftLen;f2 = halfLen * fs / fftLen;longChirp = chirp2(t, D, f1, f2)';
input = [longChirp, circshift(longChirp, 10), circshift(longChirp, 100), [1; zeros(m - 1, 1)]];
target = filter(tarB, tarA, input);
%% My implementation
boundedZero = 1;
if boundedZero
    while(1)
        b0 = randn(1, nSections);
        rt1 = roots(b0);
        if any(abs(rt1) < 1)
            break;
        end
    end
    while(1)
        b1 = randn(1, nSections);
        rt1 = roots(b1);
        if any(abs(rt1) < 1)
            break;
        end
    end
    while(1)
        b2 = randn(1, nSections);
        rt1 = roots(b2);
        if any(abs(rt1) < 1)
            break;
        end
    end
else
    b0 = randn(1, nSections);
    b1 = randn(1, nSections);
    b2 = randn(1, nSections);
end
a1 = zeros(1, nSections);
a2 = zeros(1, nSections);
for idx = 1 : nSections
    while(1)
        a = [1, randn(1, 2)];
        rt1 = roots(a);
        if any(abs(rt1) > 1)
            continue;
        else
            break;
        end
    end
    a1(idx) = a(2);
    a2(idx) = a(3);
end
% [tarB, tarA] = cheby1(6,5, 0.9);
% sos = tf2sos(tarB, tarA);
% sos = rand(size(sos));
% b0 = sos(:, 1).';
% b1 = sos(:, 2).';
% b2 = sos(:, 3).';
% a1 = sos(:, 5).';
% a2 = sos(:, 6).';
c1 = a1 + 2;
c2 = (1 + a1 + a2) ./ c1;
d1 = (2 * b0 + b1) ./ c1;
d2 = (b0 + b1 + b2) ./ (c1 .* c2);
poles = [b0(:); d1(:); d2(:); c1(:); c2(:)];
func = @(x) svfLoss(x, input, target, nSections, coefNum, parallelFilter, reverseFilter);
opt.optTol = eps * 0.01;
opt.progTol = eps * 0.01;
opt.MaxIter = 1500;
opt.MaxFunEvals = 10000;
opt.Method = 'lbfgs'; % pcg, lbfgs
weightsOpt3 = minFunc(func, poles, opt);
[f1, J1, optedsig] = func(weightsOpt3);
figure(1)
plot(target);
figure(2)
plot(optedsig)
%% Construct
d0 = weightsOpt3(nSections * 0 + 1 : nSections * 1);
d1 = weightsOpt3(nSections * 1 + 1 : nSections * 2);
d2 = weightsOpt3(nSections * 2 + 1 : nSections * 3);
c1 = weightsOpt3(nSections * 3 + 1 : nSections * 4);
c2 = weightsOpt3(nSections * 4 + 1 : nSections * 5);

a2 = c1 - 2;
a3 = (c2 .* c1) - 1 - a2;
b0 = d0;
b1 = (d1 .* c1) - 2 .* d0;
b2 = (d2 .* c1 .* c2) - b0 - b1;

sos(:, 1) = b0;sos(:, 2) = b1;sos(:, 3) = b2;sos(:, 4) = 1;sos(:, 5) = a2;sos(:, 6) = a3;
function [mse, grad, out] = svfLoss(poles, input, target, nSections, coefNum, parallelFilter, reverseFilter)
poles = ADNode(poles);
b0 = reshape(poles(coefNum * nSections * 0 + 1 : coefNum * nSections * 1), [coefNum, nSections]);
d1 = reshape(poles(coefNum * nSections * 1 + 1 : coefNum * nSections * 2), [coefNum, nSections]);
d2 = reshape(poles(coefNum * nSections * 2 + 1 : coefNum * nSections * 3), [coefNum, nSections]);
c1 = reshape(poles(coefNum * nSections * 3 + 1 : coefNum * nSections * 4), [coefNum, nSections]);
c2 = reshape(poles(coefNum * nSections * 4 + 1 : coefNum * nSections * 5), [coefNum, nSections]);
tmp2 = svfFilter(input, b0, d1, d2, c1, c2, 0, 0, parallelFilter, reverseFilter);
sub = tmp2 - target;
out = tmp2.value;
mse1 = mean(sub(:).^2);
mse = mse1.value;
grad = mse1.backprop(1);
end
function x = chirp2(t,t1,f0,f1)
beta = (f1-f0)./t1;
x = cos(2*pi * ( 0.5* beta .* (t .* t) + f0 * t));
end