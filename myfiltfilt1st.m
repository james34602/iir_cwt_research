function yout = myfiltfilt1st(b, a, x)
n = length(x);
nfact = 3; % length of edge transients
zi = (b(2) - b(1) * a(2)) / (1 + a(2));

yout = zeros(n,1);
ytemp = zeros(n+2*nfact,1);
% Single channel, data explicitly concatenated into one vector
ytemp(:) = [2*x(1)-x(nfact+1:-1:2); x; 2*x(end)-x(end-1:-1:end-nfact)];

% filter, reverse data, filter again, and reverse data again
ytemp = filter(b,a,ytemp,zi*ytemp(1),1);
ytemp = ytemp(end:-1:1);
ytemp = filter(b,a,ytemp,zi*ytemp(1),1);

% retain reversed central section of y
yout(:) = ytemp(end-nfact(1,1):-1:nfact(1,1)+1);
end