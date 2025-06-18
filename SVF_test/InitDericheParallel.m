function [Bmfwd, Bmbwd, Am, divareaSum, y2] = InitDericheParallel(s, intendedArrayLen)
if s < 0.01
    s = 0.01; 
end
% a1=1.6800;a2=-0.6803;b1=3.7350;b2=-0.2598;y1=1.7830;y2=1.7230;w1=0.6318;w2=1.9970;% 0.21 - 0.9
opt1 = [1.6800, -0.6803, 3.7350, -0.2598, 1.7830, 1.7230, 0.6318, 1.9970];
% 0.9 - 2000
opt2 = [1.68335586837891,-0.684751233906089,3.70872997040196,-0.257745603161696,1.78098256200079,1.72544489520005,0.631794849607942,1.99275458858079];
% 0.0 - 0.21
opt3 = [1.68797372365328,-0.668415608466978,3.70569106695914,-0.260419852690742,1.90158005795746,1.74193623473431,0.616888363231137,1.99534058581487];
if s < 0.21
    opt = opt3;
elseif s < 0.9
    opt = opt1;
else
    opt = opt2;
end
a1=opt(1);a2=opt(2);b1=opt(3);b2=opt(4);y1=opt(5);y2=opt(6);w1=opt(7);w2=opt(8);
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
r=close_form_residuez(bfwd, p);
Bmfwd = rpk2parf(r, p);
if ~isreal(Bmfwd)
    Bmfwd = real(Bmfwd);
end
% Anti causal part
an1 = ap1 - bp1*ap0;
an2 = ap2 - bp2*ap0;
an3 = ap3 - bp3*ap0;
an4 = - bp4*ap0;
bbwd = [an1, an2, an3, an4];
r=close_form_residuez(bbwd, p);
Bmbwd = rpk2parf(r, p);
if ~isreal(Bmbwd)
    Bmbwd = real(Bmbwd);
end
x = zeros(intendedArrayLen, 1);
x(getFFTHalfLen(intendedArrayLen)) = 1;
yp1 = filter(bfwd, a, x);
ym1 = flipud(filter(bbwd, a, [0; flipud(x)]));
ym1(1) = [];
y2 = yp1+ym1;
divareaSum = 1 / sum(y2);
end
function r = close_form_residuez(b, p)
r=[(b(1)*p(1)*p(1)*p(1) + b(2)*p(1)*p(1) + b(3)*p(1) + b(4))/((p(1) - p(2))*(p(1) - p(3))*(p(1) - p(4)))
-(b(1)*p(2)*p(2)*p(2) + b(2)*p(2)*p(2) + b(3)*p(2) + b(4))/((p(1) - p(2))*(p(2) - p(3))*(p(2) - p(4)))
(b(1)*p(3)*p(3)*p(3) + b(2)*p(3)*p(3) + b(3)*p(3) + b(4))/((p(1) - p(3))*(p(2) - p(3))*(p(3) - p(4)))
-(b(1)*p(4)*p(4)*p(4) + b(2)*p(4)*p(4) + b(3)*p(4) + b(4))/((p(1) - p(4))*(p(2) - p(4))*(p(3) - p(4)))];
end
function Bm = rpk2parf(r, p)
Bm=[r(1) + r(2),  r(3) + r(4); - p(1)*r(2) - p(2)*r(1), - p(3)*r(4) - p(4)*r(3)];
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end