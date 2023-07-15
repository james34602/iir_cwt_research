wndLen = 32;
n = wndLen + 1;
half = (n+1)/2;
x = (0:half-1)'/(n-1);
Pi = sym('pi');
w = 0.5 - 0.5*cos(2*Pi*x);
w = [w; w(end-1:-1:2)];

spectrumSymbolic = feval(symengine, 'numeric::fft', w, 'Symbolic');
spectrumSymbolic = spectrumSymbolic(1 : half);
res = subs(spectrumSymbolic, Pi, pi)