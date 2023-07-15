% %%%%%%%%%%%%%%%%%%%%%%%%% Version Dec 17, 2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input must be of the form
% f_min, f_max, sr, bpo, xNmax, op, memLimGB, w
% sr = sampling rate
% bpo = bins per octave
% xNmax decides the length of data vector "x" as a multiple of largest
% window length
% op = overlap percent
% Length of x = (xNmax) times max(Nvec)
% memLimGB = memory limit for the KK matrix (program aborts if exceeded)
% w = window to use: 1=rectwin, 2=hamming, 3=hann, etc.
%
%
% KK is a cell array convertible to a matrix using cell2mat()
% Use A = cell2mat(KK) to get the mxn transform matrix A
%
% Empty matrix output if m<n
%
% Nvec is a vector containing constant-Q window lengths calculated using
% the frequencies and the bins per octave.
%
% The KK*x matrix contains the lower freq coefficient (fewer in number)
% followed by higher freq coefficients (larger in number)
%
% f = vector of frequencies analyzed
% timeMap = Cell array containing effective time values where the analysis
% is being done depending on the window location
%
% Fixed rounding errors that existed in cqker2 and cqker3. Now we slide the
% window through "fractional" sample values and then round each value
% at the end instead of rounding off at each step before calculating the
% next one. This prevents rounding errors from accumulating progressively.
%
% May 11, 2011: Made K sparse to save memory
%
% Dec 17, 2011: Added documentation comments throughout the code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

function [KK, memreqd, Nvec, f, timeMap] = cqker4(f_min, f_max, sr, bpo, xNmax, op, memLimGB, w)

if xNmax<1
    xNmax = 2;
end

switch(w)
    case 1
        name = '@rectwin';
    case 2
        name = '@hamming';
    case 3
        name = '@hann';
    case 4
        name = '@gausswin';
    case 5
        name = '@blackman';
    case 6
        name = '@bartlett';
    otherwise
        name = '@hamming';
end

fprintf('[Running with: %d, %d, %d, %d, %1.1f, %2.1f, %d, %s]', f_min, f_max,...
    sr, bpo, xNmax, op, memLimGB, name(2:end))

if f_max > sr/2
    f_max = sr/2;
    fprintf('\n[WARNING: f_max reduced to sr/2]')
end
r = 2^(1/bpo);
Q = 1 / (r - 1);
B = floor(log(f_max/f_min) / log(r)); % number of bins
k = 0:B-1;
f = f_min*r.^k;
Nvec = round(Q*sr./f);
xlen = round(xNmax*max(Nvec));
Nlen = length(Nvec);
jp = 100-op;
rows = floor((xlen*ones(1, Nlen) - Nvec) ./ (Nvec*jp/100));

KK={0}; timeMap={0};

memreqd = 16*sum(rows)*xlen/(1024)^3/6; % GB for sparse matrix
fprintf('\n[Memory required: %d GB]', (memreqd))

if sum(rows)<xlen
    disp '[m<n, not calculating.]';
end

KK = cell(Nlen, 1);
timeMap = KK;
for i=1:Nlen
    curwinstart = [];
    edge = 1;
    while(1)
        curwinstart = [curwinstart, edge];
        edge = edge+Nvec(i)*jp/100;
        if round(edge)>xlen-Nvec(i)
            break;
        end
    end
    curwinstart = round(curwinstart);

    K = zeros(length(curwinstart), xlen);
    tmmp = zeros(1, length(curwinstart));

    e = exp(-1i*2*pi*Q*(0:Nvec(i)-1)/Nvec(i)).*...
        window(eval(name), Nvec(i))';
    e = e / norm(window(eval(name), Nvec(i))); % Normalize energy

    for j = 1:size(K,1)
        K(j, curwinstart(j):curwinstart(j)+Nvec(i)-1) = e;
        tmmp(j) = (curwinstart(j)-1+Nvec(i)/2)/sr;
    end

    KK{i} = sparse(K);
    timeMap{i} = tmmp;
end

fprintf('\n')

end