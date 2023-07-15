halfLen = 8193;
frameSize = (halfLen - 1) * 2;
fb = cwtfilterbank('SignalLength',halfLen - 1,'Wavelet','amor','VoicesPerOctave',48,'SamplingFrequency',1);
w1 = freqz(fb);
% w1(:, halfLen+1:frameSize) = conj(w1(:, halfLen-1:-1:2));
w1(:, halfLen+1:frameSize) = 0;
w1 = ifftshift(ifft(w1.'), 1);

fb = cwtfilterbank('SignalLength',halfLen - 1,'Wavelet','Morse','VoicesPerOctave',48,'SamplingFrequency',1);
w2 = freqz(fb);
w2(:, halfLen+1:frameSize) = conj(w2(:, halfLen-1:-1:2));
w2(:, halfLen+1:frameSize) = 0;
w2 = ifftshift(ifft(w2.'), 1);

fb = cwtfilterbank('SignalLength',halfLen - 1,'Wavelet','bump','VoicesPerOctave',48,'SamplingFrequency',1);
w3 = freqz(fb);
w3(:, halfLen+1:frameSize) = conj(w3(:, halfLen-1:-1:2));
w2(:, halfLen+1:frameSize) = 0;
w3 = ifftshift(ifft(w3.'), 1);