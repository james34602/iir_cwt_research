[Bmfwd, Bmbwd, Am] = InitDerichePrototype(100);
function [Bmfwd, Bmbwd, Am] = InitDerichePrototype(s)
% a1=1.6800;a2=-0.6803;b1=3.7350;b2=-0.2598;y1=1.7830;y2=1.7230;w1=0.6318;w2=1.9970;
a1=1.6800;a2=-0.679999;b1=3.71;b2=-0.2598;y1=1.7830;y2=1.7230;w1=0.6318;w2=1.9970;
syms a1 a2 b1 b2 y1 y2 w1 w2 s
% Causal part
ap0 = a1 + a2;
ap1 = exp(-y2/s)*(b2*sin(w2/s) - (a2 + 2*a1)*cos(w2/s)) + exp(-y1/s)*(b1*sin(w1/s)- (2*a2+a1)*cos(w1/s));
ap2 = 2*exp(-(y1+y2)/s)*((a1+a2)*cos(w2/s)*cos(w1/s) - cos(w2/s)*b1*sin(w1/s) - cos(w1/s)*b2*sin(w2/s)) + a2*exp(-2*y1/s) + a1*exp(-2*y2/s);
ap3 = exp(-(y2+2*y1)/s)*(b2*sin(w2/s) - a2*cos(w2/s)) + exp(-(y1 + 2*y2)/s)*(b1*sin(w1/s) - a1*cos(w1/s));
bp4 = exp(-(2*y1 + 2*y2)/s);
bp3 = - 2*cos(w1/s)*exp(-(y1+2*y2)/s) - 2*cos(w2/s)*exp(-(y2+2*y1)/s);
bp2 = 4*cos(w2/s)*cos(w1/s)*exp(-(y1+y2)/s) + exp(-2*y1/s) + exp(-2*y2/s);
bp1 = - 2*exp(-y2/s)*cos(w2/s) - 2*exp(-y1/s)*cos(w1/s);
bfwd = [ap0, ap1, ap2, ap3];
a = [1, bp1, bp2, bp3, bp4];
%%
p = [(sqrt(-exp(2*y2/s))*sin(w2/s) + exp(y2/s)*cos(w2/s))*exp(-2*y2/s)
(-sqrt(-exp(2*y2/s))*sin(w2/s) + exp(y2/s)*cos(w2/s))*exp(-2*y2/s)
(sqrt(-exp(2*y1/s))*sin(w1/s) + exp(y1/s)*cos(w1/s))*exp(-2*y1/s)
(-sqrt(-exp(2*y1/s))*sin(w1/s) + exp(y1/s)*cos(w1/s))*exp(-2*y1/s)];
Am = [1, 1; - p(1) - p(2), - p(3) - p(4); p(1)*p(2), p(3)*p(4)];
p = sym('p', [4, 1]);
r=residued2(bfwd, p); % perform partial fraction expansion to the delayed form
Bmfwd = rpk2parf(r, p); % recombine to second-order sections
% [r,p,f2] = residuez(bfwd,a);
% [Bm,Am_,FIR]=rpk2parf2(r,p,[], 1); % recombine to second-order sections
if ~isreal(Bmfwd)
    Bmfwd = real(Bmfwd);
end
% Anti causal part
an1 = ap1 - bp1*ap0;
an2 = ap2 - bp2*ap0;
an3 = ap3 - bp3*ap0;
an4 = - bp4*ap0;
bbwd = [an1, an2, an3, an4];
r=residued2(bbwd, p); % perform partial fraction expansion to the delayed form
Bmbwd = rpk2parf(r, p); % recombine to second-order sections
if ~isreal(Bmbwd)
    Bmbwd = real(Bmbwd);
end
end
function r = residued2(b, p)
r=[(b(1)*p(1)*p(1)*p(1) + b(2)*p(1)*p(1) + b(3)*p(1) + b(4))/((p(1) - p(2))*(p(1) - p(3))*(p(1) - p(4)))
-(b(1)*p(2)*p(2)*p(2) + b(2)*p(2)*p(2) + b(3)*p(2) + b(4))/((p(1) - p(2))*(p(2) - p(3))*(p(2) - p(4)))
(b(1)*p(3)*p(3)*p(3) + b(2)*p(3)*p(3) + b(3)*p(3) + b(4))/((p(1) - p(3))*(p(2) - p(3))*(p(3) - p(4)))
-(b(1)*p(4)*p(4)*p(4) + b(2)*p(4)*p(4) + b(3)*p(4) + b(4))/((p(1) - p(4))*(p(2) - p(4))*(p(3) - p(4)))];
end
function Bm = rpk2parf(r, p)
Bm=[r(1) + r(2),  r(3) + r(4); - p(1)*r(2) - p(2)*r(1), - p(3)*r(4) - p(4)*r(3)];
end