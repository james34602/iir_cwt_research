o = emaFFTSpec_m(data,fs,4096,32,6,5);
magSpec = abs(o)';
imagesc(power(magSpec, 0.1))
colormap(jet);
set(gca,'YDir','normal');
function o = emaFFTSpec_m(y,fs,windowSize,stepSize,fgt_facC,m)
freqs = 0:fs./windowSize:fs./2;
freqs(1) = 1;
fgt_facT = 1./freqs.*fgt_facC;
fgt_fac = 1-10.^(log10(0.5)./fgt_facT.*stepSize./fs);
fgt_fac = fgt_fac.';
temp = [zeros(windowSize/2,1)
    y
    zeros(windowSize/2-1,1)];
lennon = windowSize/2+1;
test = zeros(windowSize,1);
test(stepSize+1) = 1;
ang = angle(fft(test));
ang = ang(1:lennon);
angC = 0;
o = zeros(ceil(length(y)./stepSize),lennon);
window1 = zeros(windowSize,1);
len = (windowSize-stepSize.*2)./2;
window1(len+1:len+stepSize.*2) = hann(stepSize.*2, 'periodic');
c = 1;
curAng = zeros(lennon,1);
curMag = zeros(lennon,1);
ema = zeros(lennon,m);
for a = 1:stepSize:length(temp)-windowSize
    temp2 = fft(temp(a:a+windowSize-1).*window1);
    temp2 = temp2(1:lennon);
    curAng(:) = angle(temp2)+angC;
    curMag(:) = abs(temp2);
    ema(:,1) = (cos(curAng)+i*sin(curAng)).*curMag.*fgt_fac+ema(:,1).*(1-fgt_fac);
    for b = 2:m
      ema(:,b) = ema(:,b-1).*fgt_fac+ema(:,b).*(1-fgt_fac);
    end
    o(c,:) = abs(ema(:,m));
    angC = angC+ang;
    angC(angC<-2.*pi) = angC(angC<-2.*pi)+2.*pi;
    c = c+1;
end
end