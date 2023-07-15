function [f, g, h] = fitExp(corrF, cdf, raise)
f = mean( ( (corrF .^ raise) - cdf) .^ 2, 1); % Cost
if nargout > 1
    J1 = -(2*corrF.^raise.*log(corrF) .* (cdf - corrF.^raise)); % Gradient
    J1(isnan(J1)) = 0;
    g = mean(J1, 1);
end
if nargout > 2
    H1 = 4 * [corrF.^(2*raise).*log(corrF).^2; -corrF.^raise.*log(corrF).^2.*(cdf - corrF.^raise)]; % Hessian
    H1(isnan(H1)) = 0;
    h = mean(H1, 1);
end
end