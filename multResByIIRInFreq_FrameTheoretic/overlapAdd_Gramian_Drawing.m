function overlapAdd_Gramian_Drawing(N, hop)
if nargin < 1
    N = 64;
    hop = 16;
end
dftMtx = dftmtx(N);
gTime = hann(N, 'periodic').^5.5 .';
windowedDFTMatrix2 = dftMtx .* gTime;
LSsolution = ifft(windowedDFTMatrix2.');
Phi0 = LSsolution * dftMtx;
[M, N] = size(Phi0);

ckmin = 2;
ck = max(ckmin,floor(N/M)); % determines size of matrix Phi (and reconstruction error)
Nx = ck*N;
Nl = floor(Nx/hop);
if hop*(Nl-1) + N < ck*N
    disp([hop*(Nl-1) + N, ck*N])
    error('Synthesis frame matrix does''t exist');
end
tic
Phi_ = sparse(M*Nl,(Nl-1)*hop+N);
for m = 0:Nl-1
    Phi_(m*M+1:m*M+M,m*hop+1:m*hop+N) = Phi0;
end
tPhi_ = Phi_';
C = real(tPhi_*Phi_);
C2 = C(N+1:ck*N, (N+1:ck*N));
figure(1)
clf
subplot(2, 1, 1)
correctionWndHF = overlapAdd2(repmat(gTime.^2, Nl, 1 )', hop);
axis tight
title('Overlap-add illustration')
subplot(2, 1, 2)
plot(diag(C) / N);
hold on
plot(correctionWndHF)
hold off;
axis tight
legend('Frame operator Gramian diagonal component', 'Overlap-add squared')
title('Gramian diagonal vs Overlap-add')
end
function y = overlapAdd2(tmp, hop)
nframes = size(tmp, 2);
fftLen = size(tmp, 1);
xlen = fftLen + (nframes-1)*hop;
y = zeros(xlen, 1);
plt = zeros(xlen, 1);
hold on
for l = 1 : nframes
    plt(:) = 0;
    plt(1+(l-1)*hop : fftLen+(l-1)*hop) = tmp(:, l);
    plot(plt)
    y(1+(l-1)*hop : fftLen+(l-1)*hop) = y(1+(l-1)*hop : fftLen+(l-1)*hop) + tmp(:, l);
end
hold off
end