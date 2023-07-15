rng(1)
clear
m = 4;
nSigs = 3;
input = randn(m, nSigs);
dp = randn(nSigs, 1);
nSections = 2;
ltv = 1;
if ltv == 1
    coefNum = m;
else
    coefNum = 1;
end
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
poles = [input(:); b0(:); d1(:); d2(:); c1(:); c2(:)];
if ~statesNoGrad
    poles = [poles; initialValueZ1(:); initialValueZ2(:)];
end
poles = ADNode(poles);
input = reshape(poles(1 : m * nSigs), [m, nSigs]);
b0 = reshape(poles(m * nSigs + coefNum * nSections * 0 + 1 : m * nSigs + coefNum * nSections * 1), [coefNum, nSections]);
d1 = reshape(poles(m * nSigs + coefNum * nSections * 1 + 1 : m * nSigs + coefNum * nSections * 2), [coefNum, nSections]);
d2 = reshape(poles(m * nSigs + coefNum * nSections * 2 + 1 : m * nSigs + coefNum * nSections * 3), [coefNum, nSections]);
c1 = reshape(poles(m * nSigs + coefNum * nSections * 3 + 1 : m * nSigs + coefNum * nSections * 4), [coefNum, nSections]);
c2 = reshape(poles(m * nSigs + coefNum * nSections * 4 + 1 : m * nSigs + coefNum * nSections * 5), [coefNum, nSections]);
if ~statesNoGrad
    initialValueZ1 = reshape(poles(m * nSigs + coefNum * nSections * 5 + 1 : m * nSigs + coefNum * nSections * 5 + nSigs * nSections), [nSigs, nSections]);
    initialValueZ2 = reshape(poles(m * nSigs + coefNum * nSections * 5 + nSigs * nSections + 1 : m * nSigs + coefNum * nSections * 5 + nSigs * nSections * 2), [nSigs, nSections]);
else
    initialValueZ1 = 0;
    initialValueZ2 = 0;
end
tmp2 = svfFilter(input, b0, d1, d2, c1, c2, initialValueZ1, initialValueZ2, parallelFilter, reverseFilter);
sub = tmp2 * dp - target;
mse1 = sum(sub.^2);
grad3 = mse1.backprop(1);
plot(grad3);axis tight