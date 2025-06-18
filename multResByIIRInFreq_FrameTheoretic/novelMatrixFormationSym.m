function novelMatrixFormationSym(fftLen, hop)
if nargin < 1
    fftLen = 3;
    hop = 2;
end
Phi0 = randn(fftLen, fftLen);
useSym = 1;
if useSym
    Phi0 = sym('a', [fftLen, fftLen], 'real');
end
[halfLen, fftLen] = size(Phi0);
% build frame matrix Phi of Xa = Phi*x
% Phi is used to compute pseudo inverse (dual frame matrix)
ckmin = 3;
ck = max(ckmin,floor(fftLen/halfLen)); % determines size of matrix Phi (and reconstruction error)
Nx = ck*fftLen;
Nm = floor(Nx/hop);
k = ck*floor(fftLen/hop);
ic = halfLen*floor(k*halfLen/(2*halfLen)); % index to select matrix block
if hop*(Nm-1) + fftLen < ck*fftLen
    disp([hop*(Nm-1) + fftLen, ck*fftLen])
    error('Synthesis frame matrix does''t exist');
end
% tic
Phi_ = zeros(halfLen*Nm,(Nm-1)*hop+fftLen);
if useSym
    Phi_ = sym(Phi_);
end
for m = 0:Nm-1
    Phi_(m*halfLen+1:m*halfLen+halfLen,m*hop+1:m*hop+fftLen) = Phi0;
end
C_ = Phi_'*Phi_; % Gram without chopping head and tail
C2_ = C_(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen)); % Gram chopped
% toc
tic
blk = Phi0'*Phi0;
w = size(blk,2);
h = size(blk,1);
N = Nm;
vshft = hop;
hshft = hop;
nRows = vshft*(N-1) + h;
nCols = hshft*(N-1) + w;
% C = sparse(nRows,nCols);
C = zeros(nRows,nCols);
if useSym
    C = sym(C);
end
% [jdx0,idx0] = meshgrid(1:w,1:h);
for l = 2:N
    % idx = idx0 + (l-1)*vshft;
    % jdx = jdx0 + (l-1)*hshft;
    % in case you want really overlapping to add overlapping blocks
    % if l > 2
        % C = C + sparse(idx(:),jdx(:),blk(:),nRows,nCols);
        % C((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) = C((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) + blk;
    % else
        % C = sparse(idx(:),jdx(:),blk(:),nRows,nCols);
    %     C((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) = blk;
    % end
    C((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) = C((1:h) + (l-1)*vshft, (1:w) + (l-1)*hshft) + blk;
    C2 = C(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen));
end
C2 = C(fftLen+1:ck*fftLen, (fftLen+1:ck*fftLen));
dif = C2 - C2_;
tf=isequal(C2, C2_)
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function spec = ltv_1st_2ndNoSq2(dftSpec, b, a, c1, c2)
%% Gaussian windowing
tmp = zeros(size(dftSpec, 1), 1, 'like', dftSpec);
spec = ltv1Slice(dftSpec, tmp, b, a, c1, c2);
end