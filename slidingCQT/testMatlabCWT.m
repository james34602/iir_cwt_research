[x, fs] = loadSignal(7, -1);
x = x(1 : 3000);

[wt, cwtF] = cwt(x,fs,'amor','VoicesPerOctave',48);
wt = wt ./ max(abs(wt(:)));
h = pcolor(0 : (length(x) - 1), cwtF, 20 * log10(abs(wt)));
set(h,'EdgeColor', 'none')
colorbar
xrec = icwt(wt,'amor','SignalMean',mean(x),'VoicesPerOctave',48);
10*log10(sum(abs(x).^2)/sum(abs(x-xrec(:)).^2))