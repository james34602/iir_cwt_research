% mex -R2018a -g -output gridMapping gridMapping.c
% mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -DNOCHECK -output gridMapping gridMapping.c
% [data, fs] = audioread('test_audioCut.wav');
a = 0.00000000;
lb = 0.01;
ub = fs / 2;
% identity = eye(halfLen);
% transformationMatrix = zeros(halfLen);
% for ii = 1 : halfLen
%     [y, map] = gridMapping(identity(:,ii), a, lb, ub, fs);
%     x = interp1(map(:, 1), y, 0 : (halfLen - 1), 'makima');
%     transformationMatrix(:,ii) = x;
% end
% imagesc(transformationMatrix)
frameLen = 1024;
halfLen = frameLen / 2 + 1;
wnd = hann(frameLen, 'periodic');
hop = 512;
paddedSignal = [zeros(frameLen - hop, 1); data];
nframes = fix(ceil(length(paddedSignal) / hop) - frameLen / hop / 2);
frameindex = 1 + (0 : nframes - 1) * hop;
loop_i = 0;
cnt = 0;
for ii = 1 : size(data, 1)
    if cnt == hop
        loop_i = loop_i + 1;
        frame = (paddedSignal(frameindex(loop_i) : frameindex(loop_i) + frameLen - 1)) .* wnd;
        spec = fft(frame);
        spec = spec(1 : halfLen);
        [y, map] = gridMapping(spec, a, lb, ub, fs);
        [ySpec, map] = gridMapping(abs(spec), a, lb, ub, fs);
        xSpec = interp1(map(:, 1), ySpec, 0 : (halfLen - 1), 'makima');
        plot(abs(spec));
        hold on;
        plot(xSpec);
        hold off;
        axis tight
        cnt = 0;
    end
    cnt = cnt + 1;
end