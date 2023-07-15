% Wrapper for running the phase vocoding algorithm implemented in cqphvoc.m
% Inputs: filename (.wav), stretch factor, LSICQS parameters
% Output: Writes the edited file to filename_edited.wav


filename = 'Alto_Sax_08.wav'; scale = 1;
try
    [soundfile,sr0] = audioread(filename);
catch ME
    disp(['Error: ', ME.message]);
    return;
end
% Set all the parameters for the LSICQS calculation
sr = 8000;
bpo = 24;
f_min = 110;
f_max = 7500;
op = 95;
memlim = 4;
wind = 2;
g = gcd(sr0, sr);

% Design noise filter for removing spurious signals < f_min Hz
h  = fdesign.highpass(round(f_min-f_min/10), f_min, 80, 1, sr);
Hd_hpf = design(h, 'butter', 'MatchExactly', 'stopband');

% Downsample audiofile for faster processing and memory efficiency
soundfile = resample(soundfile, sr/g, sr0/g);
soundfile = filter(Hd_hpf, soundfile);
audiowrite([filename(1:end-4),'_original.wav'], soundfile, sr);

% Form the transform matrix A
[kk, ~, Nvec, f, timeMap] = cqker4(f_min, f_max, sr, bpo, 3, op, memlim, wind);
disp('[Transform matrix ready.]');
opprime = 1-scale*(1-op/100);

% Form the inversion matrix A'
[kk1p5, ~, ~, ~, timeMap1p5] = cqker4(f_min, f_max, sr, bpo, ...
    3*scale, opprime*100, memlim, wind);
kk1p5reim = sparse([real(cell2mat(kk1p5)); imag(cell2mat(kk1p5))]);
disp('[Inverse matrix ready.]');

% Set frame sizes for input and output signals
seglen_in = 3*max(Nvec);
seglen_out = size(cell2mat(kk1p5),2); % approx equal to scale*seglen_in
warning off all % Supress warnings from mldivide

% Fix indices of the left edge of input and output frames for overlap-add
curwinstart_in = [];
curwinstart_out=[];
edge_in = 1;
edge_out=1;
jump = 90;
while(1)
    curwinstart_in = [curwinstart_in, edge_in];
    edge_in = edge_in+seglen_in*jump/100;
    if round(edge_in)>length(soundfile)-seglen_in
        break;
    end
end
buildsignal = zeros(seglen_out*floor(length(soundfile)/seglen_in), 1);
while(1)
    curwinstart_out = [curwinstart_out, edge_out];
    edge_out = edge_out+seglen_out*jump/100;
    if round(edge_out)>length(buildsignal)-seglen_out
        break;
    end
end
curwinstart_in = round(curwinstart_in);
curwinstart_out = round(curwinstart_out);

% Process each frame with the LSICQS phasevocoder, (using Tukey Window)
for section=1:min(length(curwinstart_in),length(curwinstart_out))
    disp(['In section ', num2str(section), '/', num2str(length(curwinstart_in))]);
    z = soundfile(curwinstart_in(section):curwinstart_in(section)+seglen_in-1);
    z = z.*tukeywin(length(z),2*(1-jump/100));
    CQmagph = cqphvoc(z, f, sr, bpo, op, timeMap, timeMap1p5, Nvec, kk, scale);
    CQrecons1p5 = [];
    for i=1:length(CQmagph)
        if length(CQmagph{i}) >= length(timeMap1p5{i})
            CQrecons1p5 = [CQrecons1p5; CQmagph{i}(1:length(timeMap1p5{i}))];
        else
            CQrecons1p5 = [CQrecons1p5; CQmagph{i}; CQmagph{i}(end)*ones(length(timeMap1p5{i}) - length(CQmagph{i}),1)];
        end
    end
    burst = kk1p5reim\[real(CQrecons1p5); imag(CQrecons1p5)];
    % Filter out low frequency noise (< f_min)
    burstf = filter(Hd_hpf, burst);
    % Undo effect of Tukey window
    burstf = burstf./tukeywin(length(burstf), 2*(1-jump/100));
    % Clip signal to [-1,1] to remove spikes, if any
    burstf(burstf>=1.1) = 1;
    burstf(burstf<=-1.1) = -1;
    % Reapply Tukey window
    burstf = burstf.*tukeywin(length(burstf), 2*(1-jump/100));
    % Overlap-add the burst
    buildsignal(curwinstart_out(section):curwinstart_out(section)+length(burst)-1) = ...
        buildsignal(curwinstart_out(section):curwinstart_out(section)+length(burst)-1) + burstf;

end

% Filter LF noise from the final overlap-added signal
filteredsignal = filter(Hd_hpf, buildsignal);
disp('[Bass noise filtered.]');

% Write final output to wav file
audiowrite([filename(1:end-4), '_edited.wav'], filteredsignal, sr);
disp(['Edited wav file written to ',filename(1:end-4), '_edited.wav']);