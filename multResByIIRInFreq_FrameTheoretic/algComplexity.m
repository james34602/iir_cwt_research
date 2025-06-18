function algComplexity(sigLen, fftLen, cplxInput, needReconstructSignal, fs, lgScale)
rangex = [1, 16];
rangey = [0, 1e12];
hop = (1 : 1 : 256);
complexity1 = zeros(length(hop), 3);
complexity2 = zeros(length(hop), 3);
for trial = 1 : length(hop)
    [totalComplexitystft, totalComplexity2nd, totalComplexity4th] = proposedComplexity(sigLen, fftLen, hop(trial), needReconstructSignal, cplxInput, 0);
    [totalComplexitystftpruned, totalComplexity2ndpruned, totalComplexity4thpruned] = proposedComplexity(sigLen, fftLen, hop(trial), needReconstructSignal, cplxInput, 1);
    complexity1(trial, 1) = totalComplexitystft;
    complexity1(trial, 2) = totalComplexity2nd;
    complexity1(trial, 3) = totalComplexity4th;
    complexity2(trial, 1) = totalComplexitystftpruned;
    complexity2(trial, 2) = totalComplexity2ndpruned;
    complexity2(trial, 3) = totalComplexity4thpruned;
end
subbands = (64 : 128 : 1024);
complexity3 = zeros(length(subbands), 2);
for trial = 1 : length(subbands)
    [complexity3(trial, 1), complexity3(trial, 2)] = fCWTComplexity(sigLen, subbands(trial), cplxInput, fs);
end
clf
plot(hop, complexity2);
% yline(fCWTComplexity(sigLen, 256, cplxInput))
axis tight
% clf
hold on
rng(1)
labels = {"STFT", "Second order", "Fourth order"};
for trial = 1 : length(subbands)
    if lgScale == 1
        yl = plot(1, complexity3(trial, 1), 'x','LineWidth',3);
    else
        yl = plot(1, complexity3(trial, 2), 'x','LineWidth',3);
    end
    yl.Color = [rand(), rand(), rand()];
    labels = [labels, "fCWT subbands = " + string(subbands(trial))];
end
hold off
xlim(rangex)
ylim(rangey)
legend(labels)
xlabel("Hop")
ylabel("Multiplications")
end
function [totalComplexitystft, totalComplexity2nd, totalComplexity4th] = proposedComplexity(sigLen, fftLen, hop, needReconstructSignal, cplxInput, prunedFFT)
halfLen = fftLen / 2 + 1;
if needReconstructSignal == false && cplxInput == false
    prepad = round(halfLen / 8);
    pospad = round(halfLen / 16);
    subbandLen = halfLen + prepad + pospad;
else
    subbandLen = fftLen;
end
nframes = (fftLen / hop) + sigLen / hop + (fftLen / hop);
if prunedFFT == 1
    if hop > 1
        W = exp((0 : (halfLen - 1)) * hop * (2 * pi) / fftLen * 1i);
        tol = 1000*eps;
        nonZero = (halfLen * 2) - (sum(abs(abs(real(W)) - 1) < tol)+sum(abs(abs(real(W)) - 0) < tol)+sum(abs(abs(imag(W)) - 1) < tol)+sum(abs(abs(imag(W)) - 0) < tol));
        fftComplexity = nonZero + fftLen * log2(hop);
    else
        fftComplexity = fftLen * log2(hop) + halfLen;
    end
else
    fftComplexity = fftLen * log2(fftLen);
end
perFrameComplexitystft = fftLen + fftComplexity; % Windowing + FFT

totalComplexitystft = nframes * perFrameComplexitystft;

perFrameComplexity2nd = fftComplexity + ... % FFT
subbandLen * 6 + ... % LTV-IIR forward
subbandLen * 6 + ... % LTV-IIR reverse
halfLen; % Frequency dependent weight

totalComplexity2nd = nframes * perFrameComplexity2nd;

perFrameComplexity4th = fftComplexity + ... % FFT
subbandLen * 10 + ... % LTV-IIR filter 1 forward
subbandLen * 10 + ... % LTV-IIR filter 2 forward
subbandLen * 10 + ... % LTV-IIR filter 1 reverse
subbandLen * 10 + ... % LTV-IIR filter 2 reverse
halfLen; % Frequency dependent weight

totalComplexity4th = nframes * perFrameComplexity4th;
end
function [complexity1, complexity2] = fCWTComplexity(sigLen, subbandLen, cplxInput, fs)
% sigLen = 2^nextpow2(sigLen);
fCWTConvolveLogFreqComplexity = fcwtComplexity(sigLen, fs, 0.1, fs/2, subbandLen, 0);
fCWTConvolveLinFreqComplexity = fcwtComplexity(sigLen, fs, 0.1, fs/2, subbandLen, 1);
if cplxInput == false
    totalConvolveLogComplexity = sum(fCWTConvolveLogFreqComplexity * 2); % Single sided
    totalConvolveLinComplexity = sum(fCWTConvolveLinFreqComplexity * 2); % Single sided
else
    totalConvolveLogComplexity = sum(fCWTConvolveLogFreqComplexity * 4); % Double sided
    totalConvolveLinComplexity = sum(fCWTConvolveLinFreqComplexity * 4); % Double sided
end
complexity1 = ...
sigLen * log2(sigLen) + ... % FFT
totalConvolveLogComplexity + ...
subbandLen * (sigLen * log2(sigLen));
complexity2 = ...
sigLen * log2(sigLen) + ... % FFT
totalConvolveLinComplexity + ...
subbandLen * (sigLen * log2(sigLen));
end