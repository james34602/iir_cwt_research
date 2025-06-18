addpath 'C:\Users\JamesHighPerformance\Downloads\BeSD_Parametric_Std'
addpath 'C:\Users\JamesHighPerformance\Downloads\BeSD_Parametric_Std\mex'
sigmaList = [0.1, 0.45, 0.75, 0.9, linspace(1, 200, 60), linspace(205, 1000, 120), linspace(1010, 10000, 140)];
sigmaList = smooth(smooth(sigmaList));
opt = [137 / 100, 0.5];
dim = length(opt);
finalSolution = zeros(length(sigmaList), dim);
lossList = zeros(length(sigmaList), 2);
K = 5;
N = 10;
epk = 3000;
parfor ii = 1 : length(sigmaList)
    sigma = sigmaList(ii);
    intendedArrayLen = round(sigma * 10);
    if mod(intendedArrayLen, 2)
        intendedArrayLen = intendedArrayLen + 1;
    end
    if intendedArrayLen < 512
        intendedArrayLen = 512;
    end
    imp=zeros(intendedArrayLen, 1);
    centre = getFFTHalfLen(intendedArrayLen);
    imp(centre) = 1;
    lower = [0.01, 0.01];
    upper = [5, 5];
    PRNG = james_randinit;
    initialAns = zeros( K * N , dim );
    for kk = 1 : K * N
        initialAns(kk, :) = opt + (james_rand(PRNG) * (upper - lower) + lower) * 0.1;
    end
    fnc = @(x) secondOrderOut(x, sigma, intendedArrayLen);
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
function y = secondOrderOut(fp, s, intendedArrayLen)
x = zeros(intendedArrayLen, 1);
centre = getFFTHalfLen(intendedArrayLen);
x(centre) = 1;
standardGauss = gaussmf((0:1:(intendedArrayLen-1))', [s, centre-1]);
standardGauss = standardGauss / sum(standardGauss);
y3 = zeros(size(fp, 1), intendedArrayLen);
for idx = 1 : size(fp, 1)
    opt = fp(idx, :);
    [b, a] = gauss_precompute2(opt, s);
    vYSignal = filter(b, a, x);
    vYSignal = filter(b, a, vYSignal(end:-1:1));
    y2 = vYSignal(end:-1:1);
    dif = y2 / sum(y2) - standardGauss;
    y3(idx, :) = dif .^ 2;
end
y = mean(y3, 2);
end
function [b, a] = gauss_precompute2(fp, sigma)
mp = exp(-fp(1) ./ sigma);
a1 = -2 * cos(fp(2) ./ sigma) .* mp;
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