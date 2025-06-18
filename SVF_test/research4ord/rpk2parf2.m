function [Bm,Am,FIR]=rpk2parf2(r,p,k,poleSorting)

%RPK2PARF converts the partial fraction form to parallel SOS form.
%   [Bm,Am,FIR]=RPK2PARF(r,p,k) converts the partial fraction form
%   computed by RESIDUE, RESIDUEZ, or RESIDUED to parallel second-order filters.
%   The Bm and Am parameters are computed from the residues r and poles p
%   which are paired as complex conjugate pairs. The FIR part simply equals the
%   direct terms FIR=k;
%
%   The Bm and Am matrices are containing the [b0 b1]' and [1 a0 a1]'
%   coefficients for the different sections in their columns. For example,
%   Bm(:,3) gives the [b0 b1]' parameters of the third second-order
%   section. These can be used by the filter command separatelly (e.g., by
%   y=filter(Bm(:,3),Am(:,3),x), or by the PARFILT or DELPARFILT command.
%   (The same function is used for normal and delayed parallel filters.)
%
%   Also, remember that converting a transfer function to parallel second-order
%   is not possible if there is pole mulitplicity.
%
%   For details, see:
%
%   Balázs Bank and Julius O. Smith, III, "A delayed parallel filter structure
%   with an FIR part having improved numerical properties", 136th AES
%   Convention, Preprint No. 9084, Berlin, April 2014.
%
%   http://www.mit.bme.hu/~bank/delparf
%
%   C. Balazs Bank, 2014.
if poleSorting
    psorted=cplxpair(p); % arrange into complex polepairs
else
    psorted=p;
end
rsorted=zeros(size(r));
for n=1:length(p),
    [val,ind]=min(abs(p-psorted(n))); % find the nearest original pole to the sorted one
    rsorted(n)=r(ind); % copy the corresponding residue
end;
p=psorted;
r=rsorted;

for s=1:floor(length(p)/2),
    [Btmp,Atmp]=residuez(r(s*2-1:s*2),p(s*2-1:s*2),[]);
    Bm(:,s)=Btmp.';
    Am(:,s)=Atmp.';
end;
if mod(length(p),2)>0, % so there is a first order filter
    [Btmp,Atmp]=residuez(r(s*2+1),p(s*2+1),[]);
    Bm(1,s+1)=Btmp.';
    Bm(2,s+1)=0;
    Am(1:2,s+1)=Atmp.';
    Am(3,s+1)=0;
end;

FIR=k;
end