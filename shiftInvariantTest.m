rng(1)
fs = 10e3;
wndLen = 256;
x = [zeros(1, wndLen), zeros(1, wndLen * 4)];
x(:) = 0;
x(wndLen + 1) = 1;
x(:)=randn(size(x));
x = [zeros(1, wndLen), x, zeros(1, wndLen)];
wnd = hamming(wndLen, 'periodic') * 0 + 1;
wnd = hamming(wndLen, 'periodic');
hop = wndLen - wndLen / 16;

[spec1,f,t1] = stft(x,fs,Window=wnd,OverlapLength=hop,FFTLength=wndLen);
[spec2,f,t2] = stft(circshift(x, round((wndLen-hop) / 2)),fs,Window=wnd,OverlapLength=hop,FFTLength=wndLen);
pwrSpec1 = abs(spec1);
pwrSpec2 = abs(spec2);

figure(1)
plot((pwrSpec1(wndLen / 4, :)))
hold on
plot((pwrSpec2(wndLen / 4, :)));
hold off
axis tight
S = perform_stft(x,wndLen,wndLen - hop, []);
S2 = perform_stft(circshift(x, round((wndLen-hop) / 2)),wndLen,wndLen - hop, []);
% energy of the signal
e = norm(x,'fro').^2;
% energy of the coefficients
eS = norm(abs(S),'fro').^2;
disp(e-eS)
function y = perform_stft(x, fftLen,q, options)
% perform_stft - compute a local Fourier transform
%
% Forward transform:
%   MF = perform_stft(M,fftLen,q, options);
% Backward transform:
%   M  = perform_stft(MF,fftLen,q, options);
%
%   fftLen is the width of the window used to perform local computation.
%   q is the spacing betwen each window.
%
%   MF(:,i) contains the spectrum around point (i-1)*q
%
%   A typical use, for an redundancy of 2 could be fftLen=2*q+1
%
%   options.bound can be either 'per' or 'sym'
%
%   options.normalization can be set to
%       'tightframe': tight frame transform, with energy conservation.
%       'unit': unit norm basis vectors, usefull to do thresholding
%
%   If fftLen and q are vectors, then it computes a multi-window (Gabor) STFT,
%   and the output of the forward transform is a cell array.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;


if length(fftLen)>1
    % Gabor STFT
    if length(fftLen)~=length(fftLen)
        error('fftLen and q must have the same length');
    end
    % mult-window STFT
    if ~iscell(x)
        % direct transform
        for i=1:length(fftLen)
            y{i} = perform_stft(x, fftLen(i), q(i), options);
        end
    else
        n = getoptions(options, 'n', 1, 1);
        y = zeros(n,size(x{1}, 3));
        % backward transform
        for i=1:length(fftLen)
            y = y + perform_stft(x{i}, fftLen(i), q(i), options);
        end
        y = y/length(fftLen);
    end    
    return;    
end


multichannel = getoptions(options, 'multichannel', 0);

if multichannel
    options.multichannel = 0;
    % multichannel transform
    if nb_dims(x)==3
        % backward transform
        for i=1:size(x,3)
             y(:,i) = perform_stft(x(:,:,i), fftLen,q, options);
        end
    else
        for i=1:size(x,2)
             y(:,:,i) = perform_stft(x(:,i), fftLen,q, options);
        end        
    end
    return;
end

if size(x,1)==1 || size(x,2)==1
    x = x(:); dir = 1;
    n = length(x);
else
    dir = -1;
    n = getoptions(options, 'n', 1, 1);
end

bound = getoptions(options, 'bound', 'per');
transform_type = getoptions(options, 'transform_type', 'fourier');
normalization = getoptions(options, 'normalization', 'tightframe');
window_type = getoptions(options, 'window_type', 'sin');
eta = getoptions(options, 'eta', 1);

% perform sampling
X = 1:q:n+1;
p = length(X);

if mod(fftLen,2)==1
% fftLen = ceil((fftLen-1)/2)*2+1;
    w1 = (fftLen-1)/2;
    dX = (-w1:w1)';
else
    dX = (-fftLen/2+1:fftLen/2)';
end

X1 = repmat(X, [fftLen 1]) + repmat(dX, [1 p]);
switch lower(bound)
    case 'sym'
        X1(X1<1) = 1-X1(X1<1);
        X1(X1>n) = 2*n+1-X1(X1>n);
    case 'per'
        X1 = mod(X1-1,n)+1;
end
I = X1;

% build a weight function
switch lower(window_type)
    case {'sin' 'hanning'}
%        t = linspace(-pi,pi,fftLen);
%        W = cos(t(:))+1;
        W = .5 *(1 - cos( 2*pi*(0:fftLen-1)'/(fftLen-1) ));
        W = gausswin(fftLen, 30);
    case 'constant'
        W = ones(fftLen,1);
    otherwise
        error('Unkwnown winow.');
end


%% renormalize the windows
ovp = fftLen / q;
correctionWndHF = overlapAdd(repmat((W.^2)', ovp * 2, 1 )', q);
correctionWndHF = correctionWndHF(fftLen - q + 1 : fftLen * 2 - q);
weight = 1 ./ sqrt(correctionWndHF);
Weight = W .* weight;
% Weight = W * mean(weight);
%% compute the transform
if dir==1
    y = zeros(eta*fftLen,p);
    if mod(fftLen,2)==1
        m = (eta*fftLen+1)/2; w1 = (fftLen-1)/2;
        sel = m-w1:m+w1;
    else
        m = (eta*fftLen)/2+1; w1 = fftLen/2;
        sel = m-w1:m+w1-1;
    end
    y(sel,:) = x(I) .* Weight;
    % perform the transform
    y = my_transform( y, +1, transform_type );
    % renormalize if necessary
    if strcmp(normalization, 'unit')
        y = y ./ repmat( Renorm, [1 p] );
    end
else
    if strcmp(normalization, 'unit')
        x = x .* repmat( Renorm, [1 p] );
    end
    x = my_transform( x, -1, transform_type );
    x = real( x.*Weight );
    y = zeros(n,1);
    for i=1:p
        y(I(:,i)) = y(I(:,i)) + x(:,i);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = my_transform(x,dir,transform_type)

% my_transform - perform either FFT or DCT with energy conservation.
%   Works on array of size (fftLen,fftLen,a,b) on the 2 first dimensions.

fftLen = size(x,1);
if strcmp(transform_type, 'fourier')
    % normalize for energy conservation
    if dir==1
        y = fft(x / sqrt(fftLen));
    else
        y = ifft( x*sqrt(fftLen) );
    end
elseif strcmp(transform_type, 'dct')
    for i=1:size(x,2)
        y(:,i) = perform_dct_transform(x(:,i),dir);
    end
else
    error('Unknown transform');
end
end
function v = getoptions(options, name, v, mendatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mendatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mendatory
    error(['You have to provide options.' name '.']);
end
end
function y = overlapAdd(tmp, hop)
nframes = size(tmp, 2);
fftLen = size(tmp, 1);
xlen = fftLen + (nframes-1)*hop;
y = zeros(xlen, 1);
for l = 1 : nframes
    y(1+(l-1)*hop : fftLen+(l-1)*hop) = y(1+(l-1)*hop : fftLen+(l-1)*hop) + tmp(:, l);
end
end