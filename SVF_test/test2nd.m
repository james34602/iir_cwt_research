for idx = 1 : length(sigmaList)
    intendedArrayLen = round(sigmaList(idx) * 30);
    if mod(intendedArrayLen, 2)
        intendedArrayLen = intendedArrayLen + 1;
    end
    if intendedArrayLen < 512
        intendedArrayLen = 512;
    end
    imp=zeros(intendedArrayLen, 1);
    centre = getFFTHalfLen(intendedArrayLen);
    imp(centre) = 1;
    x = zeros(intendedArrayLen, 1);
    centre = getFFTHalfLen(intendedArrayLen);
    x(centre) = 1;
    [b1, a1] = gauss_precompute2(sigmaList, finalSolution, sigmaList(idx));
    [b2, a2] = gauss_precompute3(sigmaList(idx));
    vYSignal = filter(b1, a1, x);
    vYSignal = filter(b1, a1, vYSignal(end:-1:1));
    y2 = vYSignal(end:-1:1);
    y2 = y2 / sum(y2);
    vYSignal = filter(b2, a2, x);
    vYSignal = filter(b2, a2, vYSignal(end:-1:1));
    y3 = vYSignal(end:-1:1);
    y3 = y3 / sum(y3);
    standardGauss = gaussmf((0:1:(intendedArrayLen-1))', [sigmaList(idx), centre-1]);
    standardGauss = standardGauss / sum(standardGauss);
    plot(standardGauss);
    hold on
    plot(y2);
    plot(y3);
    hold off
    axis tight
    disp([sigmaList(idx), norm(standardGauss - y2), norm(standardGauss - y3)])
end
function [b, a] = gauss_precompute2(sigmaList, finalSolution, sigma)
mp = exp(-interp1(sigmaList, finalSolution(:, 1), sigma) ./ sigma);
a1 = -2 * cos(interp1(sigmaList, finalSolution(:, 2), sigma) ./ sigma) .* mp;
%% Transfer function
a = [ones(size(sigma, 1), 1), a1, mp .* mp];
b = sum(a, 2);
end
function [b, a] = gauss_precompute3(sigma)
fnc1 = @(x) (1.1089588389083*x.^2 + -0.5103623243047*x + 0.31145287849522) ./ (x.^2 + -0.45770515586342*x + 0.56222236387365);
fnc2 = @(x) (0.78383661511709*x.^2 + 0.33298385757683*x + -0.30543389899193) ./ (x^2 + 0.4240173507895*x + -0.25174520144833);
fnc3 = @(x) (1.0865112833776*x.^2 + -3.6336391428065*x + 6.9734198996186) ./ (x.^2 + -3.2444036114946*x + 6.296939863196);
fnc4 = @(x) (0.79438895781586*x.^2 + -2.6512802647202*x + 5.1994276981943) ./ (x.^2 + -3.4071588847409*x + 6.6292426237822);
if sigma <= 1.8
    v1 = 0.88;
    v2 = 0.71;
elseif sigma <= 25.3
    v1 = fnc1(sigma);
    v2 = fnc2(sigma);
elseif sigma <= 51
    v1 = fnc3((sigma - 38.4) / 7.5);
    v2 = fnc4((sigma - 38.4) / 7.5);
end
if sigma > 51
    v1 = 1.07156;
    v2 = 0.802;
end
mp = exp(-v1 ./ sigma);
a1 = -2 * cos(v2 ./ sigma) .* mp;
%% Transfer function
a = [ones(size(sigma, 1), 1), a1, mp .* mp];
b = sum(a, 2);
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function y = gaussmf(x, params)
sig = params(1);
c = params(2);
y_val = @(x_val) exp((-(x_val - c)^2)/(2 * sig^2));
y = arrayfun(y_val, x);
end