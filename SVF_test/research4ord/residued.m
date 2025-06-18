function [r, p, f] = residued(b, a)
NUM = b(:)';
DEN = a(:)';
nb = length(NUM);
na = length(DEN);
f = [];
if na<=nb
  f = filter(NUM,DEN,[1 zeros(1,nb-na)]);
  NUM = NUM - conv(DEN,f);
  NUM = NUM(nb-na+2:end);
end;
[r,p,f2] = residuez(NUM,DEN);
