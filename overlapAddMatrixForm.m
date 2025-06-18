function [coeff, f] = precomuteLSExperiment(N, L, fs, oct, order, reqSynthesisWnd, sparCutOff, zp)
rng(1)
if nargin < 1
    fs = 48000;
    N = 64;
    L = 8;
    oct = 20;
    order = 2;
    reqSynthesisWnd = 1;
    sparCutOff = 5e-8;
    zp = 1;
end
M = N;
tol = 1e-16;
% build frame matrix Phi of Xa = Phi*x
% Phi is used to compute pseudo inverse (dual frame matrix)
ckmin = 6;
% ckmin = 2;
ck = max(ckmin,floor(N/M)); % determines size of matrix Phi (and reconstruction error)
Nx = ck*N;
Nl = floor(Nx/L);
if L*(Nl-1) + N < ck*N
    disp([L*(Nl-1) + N, ck*N])
    error('Synthesis frame matrix does''t exist');
end
% tic
S = zeros(M*Nl,(Nl-1)*L+N);
for m = 0:Nl-1
    i = m*M + (1 : M);
    j = m*L + (1 : N);
    S(i, j) = diag(hann(N, 'periodic'));
end
ovpMtxForm = ones(1, M*Nl) * S;
ovpAlg = overlapAdd2(repmat(hann(N, 'periodic')', Nl, 1 )', L);
end
function y = overlapAdd2(tmp, hop)
nframes = size(tmp, 2);
fftLen = size(tmp, 1);
xlen = fftLen + (nframes-1)*hop;
y = zeros(xlen, 1);
for l = 1 : nframes
    y(1+(l-1)*hop : fftLen+(l-1)*hop) = y(1+(l-1)*hop : fftLen+(l-1)*hop) + tmp(:, l);
end
end