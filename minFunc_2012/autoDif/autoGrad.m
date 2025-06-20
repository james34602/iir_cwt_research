function [f,g] = autoGrad(x,type,funObj,varargin)
% [f,g] = autoGrad(x,useComplex,funObj,varargin)
%
% Numerically compute gradient of objective function from function values
%
% type =
%     1 - forward-differencing (p+1 evaluations)
%     2 - central-differencing (more accurate, but requires 2p evaluations)
%     3 - complex-step derivative (most accurate and only requires p evaluations, but only works for certain objectives)

p = length(x);

if type == 1 % Use Finite Differencing
	f = funObj(x,varargin{:});
	mu = 2*sqrt(1e-12)*(1+norm(x));
	diff = zeros(p,1);
	for j = 1:p
		e_j = zeros(p,1);
		e_j(j) = 1;
		diff(j,1) = funObj(x + mu*e_j,varargin{:});
	end
	g = (diff-f)/mu;
elseif type == 3 % Use Complex Differentials
	mu = 1e-150;
	diff = zeros(p,1);
	for j = 1:p
		e_j = zeros(p,1);
		e_j(j) = 1;
		diff(j,1) = funObj(x + mu*i*e_j,varargin{:});
	end
	
	f = mean(real(diff));
	g = imag(diff)/mu;
else % Use Central Differencing
	mu = 2*sqrt(1e-12)*(1+norm(x));
	diff1 = zeros(p,1);
	diff2 = zeros(p,1);
	parfor j = 1:p
		e_j = zeros(p,1);
		e_j(j) = 1;
		diff1(j,1) = funObj(x + mu*e_j,varargin{:});
		diff2(j,1) = funObj(x - mu*e_j,varargin{:});
	end
	f = mean([diff1;diff2]);
	g = (diff1 - diff2)/(2*mu);
end

if 0 % DEBUG CODE
	[fReal gReal] = funObj(x,varargin{:});
	[fReal f]
	[gReal g]
	diff
	pause;
end