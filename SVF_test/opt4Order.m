addpath 'C:\Users\JamesHighPerformance\Downloads\BeSD_Parametric_Std'
addpath 'C:\Users\JamesHighPerformance\Downloads\BeSD_Parametric_Std\mex'
sigmaList = linspace(0.21, 0.9, 100);
opt = [1.6800, -0.6803, 3.7350, -0.2598, 1.7830, 1.7230, 0.6318, 1.9970];
dim = length(opt);
finalSolution = zeros(length(sigmaList), dim);
lossList = zeros(length(sigmaList), 2);
K = 5;
N = 10;
epk = 3000;
parfor ii = 1 : length(sigmaList)
    sigma = sigmaList(ii);
    intendedArrayLen = 512;
    imp=zeros(intendedArrayLen, 1);
    centre = getFFTHalfLen(intendedArrayLen);
    imp(centre) = 1;
    wt = 0.1;
    lower = opt;
    lower(opt > 0) = lower(opt > 0) - opt(opt > 0) * wt;
    lower(opt < 0) = lower(opt < 0) + opt(opt < 0) * wt;
    upper = opt;
    upper(opt > 0) = upper(opt > 0) + opt(opt > 0) * wt;
    upper(opt < 0) = upper(opt < 0) - opt(opt < 0) * wt;
    PRNG = james_randinit;
    initialAns = zeros( K * N , dim );
    for kk = 1 : K * N
        initialAns(kk, :) = opt + (james_rand(PRNG) * (upper - lower) + lower) * 0.00001;
    end
    fnc = @(x) DericheOut(x, sigma, intendedArrayLen);
    loss = fnc(opt);
    [gmin0, gbest0] = algo_besd2(fnc, K, N, dim, lower, upper, epk, initialAns);
    if gmin0 > loss
        finalSolution(ii, :) = opt;
    else
        finalSolution(ii, :) = gbest0;
    end
    lossList(ii, :) = [loss, gmin0];
end
function y = gaussmf(x, params)
sig = params(1);
c = params(2);
y_val = @(x_val) exp((-(x_val - c)^2)/(2 * sig^2));
y = arrayfun(y_val, x);
end
function y = DericheOut(fp, s, intendedArrayLen)
x = zeros(intendedArrayLen, 1);
centre = getFFTHalfLen(intendedArrayLen);
x(centre) = 1;
standardGauss = gaussmf((0:1:(intendedArrayLen-1))', [s, centre-1]);
standardGauss = standardGauss / sum(standardGauss);
y3 = zeros(size(fp, 1), intendedArrayLen);
for idx = 1 : size(fp, 1)
    opt = fp(idx, :);
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
    % Anti causal part
    an1 = ap1 - bp1*ap0;
    an2 = ap2 - bp2*ap0;
    an3 = ap3 - bp3*ap0;
    an4 = - bp4*ap0;
    bbwd = [an1, an2, an3, an4];
    yp1 = filter(bfwd, a, x);
    ym1 = flipud(filter(bbwd, a, [0; flipud(x)]));
    ym1(1) = [];
    y2 = yp1+ym1;
    dif = y2 / sum(y2) - standardGauss;
    y3(idx, :) = dif .^ 2;
end
y = mean(y3, 2) * 0.5 + max(y3, [], 2);
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end
function [gmin, gbest] = algo_besd2(fnc, K, N, D, low, up, maxcycle, S)
fitS = fnc(S); % LINE 3
[gmin, ind] = min( fitS ); % LINE 4
gbest = S( ind, : ); % LINE 4
% Iterative search phase
for epk = 1 : maxcycle % LINE 5
    % Generation of j0, j1, and j2
    while 1
        j1 = randperm(K * N);
        j2 = randperm(K * N);
        if sum(j1 == j2) == 0
            break;
        end
    end %  LINE 7
    j1 = j1(1 : N);
    j2 = j2(1 : N); % LINE 8
    j0 = j1; %  LINE 9
    
    % Setting sub-pattern matrix P and fitP
    P = S(j1, :);
    fitP = fitS(j1); % LINE 11
    
    % Generation of objective Vectors ; dv1
    dv1 = S(j2, :); % LINE 13
    
    % Top-N-Best pattern vectors
    [~ , index] = sort(fitS, 'ascend');
    H = S(index, :); % LINE 16
    tbest = P; % pre-memory
    for i = 1 : N
        tbest(i, :) = H( ceil( rand^3 * K * N ) , :); % LINE 16
    end
    
    % Generation of Bezier Mutation Vectors ; dv2
    while 1%  LINES 18 - 19
        j1 = randperm(N);
        j2 = randperm(N);
        if sum(j1 == 1:N) == 0 && sum(j2 == 1:N) == 0 && sum(j1 == j2) == 0
            break;
        end
    end
    dv2 = P; % pre-memory
    for i = 1 : N
        X = [P(i, :); P(j1(i), :); P(j2(i), :); tbest(i, :)];
        X = X(randperm(4), :);
%         B = bernsteinMatrix(3, rand);
        B = bernsteinMatrix3(rand);
        dv2(i , :) = B * X;
    end
    
    % Evolutionary Step Size
    a = randn(N, 1);
    b = 1 + rand(N , 1);
    c = randn(N, 1) .^ randi(7, N, 1);
    c(abs(c) > 11) = sign(c(abs(c) > 11)) * 11;
    F = a .* b .^ c; %  LINE 21
    % F = random('stable' , 1 , 0 , 1 , 0 , N , 1 ); % Users can use this line instead of line 59 for CAUCHY-SCALING OF BeSD
    % F = random('stable', 0.50 , 1 , 1 , 0 , N , 1 ); % Users can use this line instead of line 59 for LEVY-SCALING OF BeSD
    
    
    % Generation of Crossover Control Matrix ; map
    [map1 , map2] = genmap(N , D); % LINES 23 - 26
    if rand < rand % LINE 27
        map = map1;
    else
        map = map2;
    end
    
    % Generation of Trial Pattern Vectors ; T
    w1 = randn(N, 1);
    w2 = randn(N, 1);  %  LINE 29
    T = P + map .* F .* (w1 .* (dv2 - P) + w2 .* (dv1 - P)); % LINE 30
    
    % Boundary Control Mechanism ; % LINE 32
    p1 = rand;
    if p1 > 0.99
        for i = 1 : N
            for index = 1 : D
                if T(i, index) < low(index) || T(i, index) > up(index)
                    T(i, index) = rand * (up(index) - low(index)) + low(index);
                end
            end
        end
    else
        for i = 1 : N
            for index = 1 : D
                if T(i, index) < low(index)
                    T(i, index) = low(index);
                end
                if T(i, index) > up(index)
                    T(i, index) = up(index);
                end
            end
        end
    end
    
    % Update the sub-pattern matrix, P and fitP  ;
    fitT = fnc(T); % LINES 34 - 35
    ind = fitT < fitP;
    P(ind , :) = T(ind , :);
    fitP(ind)  = fitT(ind);
    
    %  Update the Global Solution, gbest
    [BestVal, index ] = min(fitP); % LINE 37
    BestP = P(index , :);
    if BestVal < gmin % LINE 38
        gmin = BestVal;
        gbest = BestP;
    end
    
    % Update the pattern matrix, S and fitS
    S(j0 , :) = P;
    fitS(j0) = fitP; %  LINE 40
end
end
% Generation of Crossover Control Matrix; map  ;  LINES 23 - 27
function [map1, map2] = genmap(N , D)
map1 = zeros(N, D); % LINE 23
map2 = zeros(N, D); % LINE 23
for j = 1 : N
    h = randperm(D);
    w = rand .^ randi(7); % LINE 25
    map1(j , h( 1 : ceil(w * D) ) ) = 1; % LINE 25
    h = randperm(D);
    w = 1 - rand .^ randi(7); % LINE 26
    map2(j , h( 1 : ceil(w * D) ) ) = 1; % LINE 26
end
end