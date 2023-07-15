function o = emaFFTSpec(y,fs,windowSize,stepSize,fgt_facC)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
fgt_facT = 1./freqs.*fgt_facC;
fgt_fac = 1-10.^(log10(0.5)./fgt_facT./fs);
fgt_fac = fgt_fac.';
temp = [zeros(windowSize/2,1)
    y
    zeros(windowSize/2-1,1)];
lennon = windowSize/2+1;
test = zeros(windowSize,1);
test(stepSize+1) = 1;
ang = angle(fft(test));
ang = ang(1:lennon);
o = zeros(ceil(length(y)./stepSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize-stepSize.*2)./2;
window1(len+1:len+stepSize.*2) = hann(stepSize.*2, 'periodic');
freqs = 0:fs/windowSize:fs/2;
freqs2 = 0:fs/(stepSize.*2):fs/2;
c = 1;
curAng = zeros(lennon,1);
curMag = zeros(lennon,1);
for a = 1:stepSize:length(temp)-windowSize
    temp2 = fft(temp(a:a+windowSize-1).*window1);
    temp2 = temp2(1:lennon);
    curAng = curAng-ang;
    temp3 = (cos(curAng)+i*sin(curAng)).*curMag.*(1-fgt_fac)+temp2.*fgt_fac;
    curMag = abs(temp3);
    curAng = angle(temp3);
    o(c,:) = curMag;
    c = c+1;
end