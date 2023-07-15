o = emaFFTSpec_m2(data,fs,4096,32,6,3);
magSpec = abs(o)';
imagesc(power(magSpec, 0.1))
colormap(jet);
set(gca,'YDir','normal');
function o = emaFFTSpec_m2(y,fs,windowSize,stepSize,fgt_facC,m)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
fgt_facT = 1./freqs.*fgt_facC;
fgt_fac = 1-10.^(log10(0.5)./fgt_facT.*stepSize./fs);
fgt_fac = fgt_fac.';
temp = [zeros(windowSize/2,1)
    y
    zeros(windowSize/2-1,1)];
lennon = windowSize/2+1;
phaseShifter = exp(-1i * stepSize * (2*pi*(0:lennon-1)/windowSize)');
angC = phaseShifter;
o = zeros(ceil(length(y)./stepSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize-stepSize.*2)./2;
window1(len+1:len+stepSize.*2) = hann(stepSize.*2, 'periodic');
c = 1;
ema = zeros(lennon,m);
for a = 1:stepSize:length(temp)-windowSize
    temp2 = fft(temp(a:a+windowSize-1) .* window1);
    temp2 = temp2(1:lennon);
    ema(:,1) = temp2 .* angC .* fgt_fac + ema(:,1) .* (1 - fgt_fac);
    for b = 2 : m
        ema(:, b) = ema(:, b - 1) .* fgt_fac + ema(:, b) .* (1 - fgt_fac);
    end
    o(c,:) = ema(:, m);
    angC = angC .* phaseShifter;
    c = c+1;
end
end
function o = idealHann(len)
for a = 1:len*4
    o(a) = cos(a/len.*pi);
end
o = o-min(o);
o = o./max(o);
o = o(len+1:2:len*3);
o = o.';
end