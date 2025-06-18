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
elseif testSig == 10
    fs = 500;
    %intermediate wave packet containing a zero signal
    t0 = (1/fs:1/fs:0.5)';
    p0 = 0*t0;

    %1st wave packet
    t1 = (1/fs:1/fs:5)';
    p1_final = vco(cos(2*pi*t1),[100 110],fs) + ...
        vco(cos(2*pi*t1),[20 22],fs) + ...
        vco(cos(2*pi*t1),[5 6],fs);
    L1 = length(p1_final)/10;
    p1_ham = hamming(L1*2);
    p1_window = [ones(L1*9,1); p1_ham(L1:(L1*2-1))];
    p1 = p1_final.*p1_window;

    %2nd wave packet
    t2 = (1/fs:1/fs:5)';
    p2_base = chirp(t2,100,5,50);
    p2_final = p2_base + chirp(t2,5,5,50);

    L2 = length(p2_final)/10;
    p2_ham = hamming(L2*2);
    p2_window = [p2_ham(1:L2); ones(L2*8,1); p2_ham(L2:(L2*2-1))];
    p2 = p2_final.*p2_window;

    %3rd wave packet
    t3 = (1/fs:1/fs:10)';
    p3_final = cos(2*pi*t3*(fs/1000)) + cos(2*pi*t3*(fs/500)) + cos(2*pi*t3*(fs/250));

    L3 = length(p3_final)/10;
    p3_ham = hamming(L3*2);
    p3_window = [p3_ham(1:L3); ones(L3*9,1)];
    p3 = p3_final.*p3_window;

    %Join all wave packets
    y = [p1; p0; p2; p0 ; p3];
else
    [y, fs] = audioread('example1.wav');
end
end