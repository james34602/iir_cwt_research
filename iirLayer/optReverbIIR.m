% rng(1)
% clear
m = 800;
ltv = 0;
if ltv == 1
    coefNum = m;
else
    coefNum = 1;
end
parallelFilter = 1;
reverseFilter = 0;
imp = audioread('variBRIR 0 0_48k.wav');
imp = imp(1:200, 1);
[Bm, Am, modresp] = fitImpToParallel(imp, 30, 0.58, length(imp));
nSections = size(Am, 2);
b0 = Bm(1, :);
b1 = Bm(2, :);
b2 = zeros(1, nSections);
a1 = Am(2, :);
a2 = Am(3, :);
% input = randn(m, nSigs);
fs=48000;
halfLen = 513;fftLen=1024;
cpLen=m;D = cpLen / fs;t = (0:cpLen-1) / fs;f1 = 1 * fs / fftLen;f2 = halfLen * fs / fftLen;longChirp = chirp2(t, D, f1, f2)';
input = [longChirp, circshift(longChirp, 10), circshift(longChirp, 100), [1; zeros(m - 1, 1)]];
target = filter(imp, 1, input);
%% My implementation
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
plot(target(:, 4));
hold on
% figure(2)
plot(optedsig(:, 4))
sum(abs(optedsig(:, 4)-target(:, 4)))
hold off
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
function [Bm, Am, modresp] = fitImpToParallel(impresp, ORDER, lambda, mdlLen)
p=wppoles(impresp(1:mdlLen),ORDER,lambda);
[Bm, Am] = parfiltdes(impresp,p);
imp=zeros(mdlLen,1);
imp(1)=1;
modresp=parfilt(Bm,Am,imp);
end
function [p]=wppoles(impresp,PNUM,lambda,WPLENGTH)
if nargin<4
    WPLENGTH=10000;
end
impresp_wp = warp_impres(impresp,lambda,WPLENGTH);
pwp=roots(prony2(impresp_wp,PNUM));
p=(pwp+lambda)./(1+lambda.*pwp);
end
function a = prony2(h, M)
hCol = h(:);
K = length(hCol)-1;
if K <= M
    hColPadded = [hCol; zeros(M-K+1,1)];
    K = M+1;
else
    hColPadded = hCol;
end
c = hColPadded(1,1);
if c == 0
    c = 1;
end
hColNorm = hColPadded/c;
Hin = toeplitz(hColNorm,[hColNorm(1,1), zeros(1,K)]);
if (K > M)
    H = Hin(:,1:M+1);
else
    H = Hin;
end
h1 = H((M+2):(K+1),1);
H2 = H((M+2):(K+1),2:(M+1));
a = [1; -H2\h1].';
end
function res = warp_impres(sig,lambda,n)
len=length(sig);
temp = [1 zeros(1,n)];
bw = [lambda 1]';
aw = [1 lambda]';
res = sig(1)*temp;
for i=2:len
	temp = filter(bw,aw,temp);
	res = res+sig(i)*temp;
end
end
function [Bm,Am]=parfiltdes(impresp,p)
p=p(abs(p)>0);
for k=1:length(p)
   if abs(p(k))>1			
      p(k)=1/conj(p(k));
   end
end
p=cplxpair(p);
pnum=length(p);
ppnum=2*floor(pnum/2);
ODD=0;
if pnum>ppnum
    ODD=1;
end
L=length(impresp);
imp=zeros(L,1);
imp(1)=1;
for k=1:2:ppnum-1
    resp=filter(1,poly(p(k:k+1)),imp);
    M(:,k)=resp;
    M(:,k+1)=[0 ;resp(1:L-1)];
end
if ODD
    resp=filter(1,poly(p(pnum)),imp);
    M(:,pnum)=resp;
end
A=M'*M;
b=M'*impresp(:);
par=A\b;
for k=1:ppnum/2
    Am(:,k)=poly(p(2*k-1:2*k)).';
    Bm(:,k)=par(2*k-1:2*k);
end
if ODD
    Am(:,k+1)=[poly(p(pnum)).'; 0];
    Bm(:,k+1)=[par(pnum); 0];
end
end
function y=parfilt(Bm,Am,x)
y=filter(Bm(:,1),Am(:,1),x);
s=size(Am);
for k=2:s(2)
    y = y + filter(Bm(:,k),Am(:,k),x);
end
end