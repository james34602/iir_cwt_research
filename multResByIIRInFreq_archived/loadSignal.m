function [y, fs] = loadSignal(testSig, frameSize)
if testSig == 0
    fs = 48000;
    y = zeros(frameSize, 1);
    y(1) = 1;
elseif testSig == 1
    fs = 48000; % Sampling frequency (samples per second)
    dt = 1/fs; % seconds per sample
    StopTime = 0.2; % seconds
    t = (0:dt:StopTime)'; % seconds
    F = 2000; % Sine wave frequency (hertz)
    y = sin(2*pi*F*t);
    F = 5000; % Sine wave frequency (hertz)
    y = y + sin(2*pi*F*t);
elseif testSig == 2
    fs = 48000;
    t = 0:1/fs:1-1/fs;
    dt = seconds(t(2)-t(1));
    x1 = chirp(t,100,t(end),20800,'quadratic');
    x2 = exp(2j*pi*3800*cos(2*pi*2*t));
    y = real((x1 + x2 + randn(size(t))/10))';
elseif testSig == 3
    fs = 48000;
    t = 0:1/fs:1;
    y = chirp(t,100,1,12800,'quadratic')' + chirp(t,13500,1,50)';
    y(30000) = 30;
elseif testSig == 4
    rng(1)
    load('../signal1')
    fs = 48000;
    y = [s];
elseif testSig == 5
    fs = 48000;
    y = randn(fs * 6, 1);
    y = y - mean(y);
    y = y ./ max(abs(y));
elseif testSig == 6
    fs = 48000;
    load batsignal;
    y = batsignal;
    y = y - mean(y);
elseif testSig == 8
    N  = 1024;
    t  = (0:N-1)/sqrt(N);
    a = exp(-8*(t/sqrt(N)-0.5).^2);
    s1 = a.*sin(2*pi*(250*t/sqrt(N)+50*(t/sqrt(N)).^3));
    s2 = a.*sin(2*pi*(130*t/sqrt(N)+100*(t/sqrt(N)).^2));
    s3 = a.*sin(2*pi*(90*t/sqrt(N)+0.2*cos(3*pi*(t/sqrt(N)))));
    y = (s1+s2+s3)';
    fs = 48000;
elseif testSig == 9
    fs = 48000;
    s = load('demoSignal/signal1.mat');
    y = s.sig';
else
    [y, fs] = audioread('../example1.wav');
end
end