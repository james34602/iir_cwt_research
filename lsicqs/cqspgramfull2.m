%%%%%%%%%%%%%%%%%%%%%%%%%%% Version Dec 17, 2010 %%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates interpolated constant Q transform by using pertinent rows from
% the spectrogram() and the correct window length as given by the constant Q
% window length calculations
%
% OUTPUTS:
% C = "triangular" constant Q spectrogram cell array -- cannot be surf()ed
% (C contains complex values only)
% timeMap = corresponding time values for C
% CQ = interpolated constant Q spectrogram amenable to displaying using surf() 
% (Interpolation is done using griddatan, Matlab builtin
%
% INPUTS:
% filename = string for .wav file
% starttime, stoptime = where the .wav file should be clipped
% bpo = bins per octave
% decifact = decimation factor, use for saving memory and faster processing
% interpmethod = string should be 'linear' or 'nearest'
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [CQ, C, timeMap] = cqspgramfull2(filename, f_min, f_max, starttime, stoptime, ...
    overlappercent, bpo, decifact, interpmethod)

[x, sr] = wavread(filename);
x = decimate(x, round(decifact)); sr = sr/round(decifact);

if length(x)/sr < 12
    disp('Warning: Segment too small. Duplicating.');
    x = repmat(x, floor(12*sr/length(x)), 1);
end


if(sr*stoptime > length(x))
    stoptime = floor(length(x)/sr);
    x = x(sr*starttime + 1 : end);
    disp('Warning: stoptime exceeds audio length. Taking full length.');
else
    x = x(sr*starttime + 1 : sr*stoptime + 1);
end


% Calculate constants
r = 2^(1/bpo);
Q = 1 / (r - 1);
B = floor(log(f_max/f_min) / log(r)); % number of bins
%Nmax = round(Q*sr/f_min); % maximum Nkcq

k = 0:B-1;
f = f_min*r.^k;

Nvec = round(Q*sr./f);  %Goes from N_{max} to N_{min}
%Nvec = fliplr(Nvec);  %Now Nvec goes from N_{min} to N_{max}
Nlen = length(Nvec);
timeMap = cell(Nlen,1);
noverlap = round(Nvec*overlappercent/100);

C = cell(Nlen, 1);
h = waitbar(0, 'Calculating spectrograms...');
for i=1:Nlen
    tic;
    waitbar(i/Nlen);
    
    % Calculate the whole spectrogram and retain the Q^th row only
    % See Eq (3) in LSICQS Paper
    [C{i,1}, ~, timeMap{i,1}, ~] = spectrogram(x, gausswin(Nvec(i)),...
        noverlap(i), Nvec(i), sr);
    C{i,1} = (C{i,1}(round(Q),:));
    timeMap{i,1} = timeMap{i,1}';
    if i==1
        time=toc;
        fprintf('Estimated time: %d seconds.', Nlen*time);
    end
end
close(h);

% Pick out only pertinent rows from each col of CQ for interpolation
fdata = [];
tdata = [];
sdata = [];

fprintf('\nShaping data for interpolation.\n');
for i=1:Nlen
    %idx1 = find(freqMap{i,1} <= f(length(Nvec)+1 - i), 1, 'first');
    %idx2 = find(freqMap{i,1} >= f(length(Nvec)+1 - i), 1, 'last');
    %if(~isempty(idx1))
        %disp('a');
        fdata = [fdata; f(i)*ones(length(timeMap{i,1})-1,1)];
        %disp('b');
        tdata = [tdata; timeMap{i,1}(1:end-1)];
        %disp('c');
        %size(CC{i,1}(idx1,:))
        sdata = [sdata, abs(C{i,1}(1, 1:length(timeMap{i,1})-1))];
    %end    
end
sdata = sdata';

% Create CQ matrix using ``pertinent'' data and griddatan() linear
if strcmp(interpmethod,'linear') || strcmp(interpmethod, 'l')
    disp('Interpolating data for final display. This may take a few minutes...');
    [tt, ff] = meshgrid(0:0.01:(stoptime-starttime-2),interp(f,4));
    %CQl = zeros(size(tt, 1), size(tt,2));
    numPoints = size(tt,1)*size(tt,2); 
    ttt = reshape(tt, numPoints, 1);
    fff = reshape(ff, numPoints, 1);
    CQ = griddatan([tdata, fdata], sdata, [ttt, fff], 'linear');
    CQ = reshape(CQ, size(tt, 1), size(tt,2));
    figure; surf((CQ).^0.4,'edgecolor','none'); axis tight; view(0,90)
    title('Linear interpolation with griddatan()');
    axis tight; view(0,90);

% Create CQ matrix using ``pertinent'' data and griddatan() nearest
elseif strcmp(interpmethod,'nearest') || strcmp(interpmethod, 'n')
    disp('Interpolating data for final display. This may take a few minutes...');
    [tt, ff] = meshgrid(0:0.01:(stoptime-starttime-2),interp(f,2));
    %CQn = zeros(size(tt, 1), size(tt,2));
    numPoints = size(tt,1)*size(tt,2); 
    ttt = reshape(tt, numPoints, 1);
    fff = reshape(ff, numPoints, 1);
    CQ = griddatan([tdata, fdata], sdata, [ttt, fff], 'nearest');  
    CQ = reshape(CQ, size(tt, 1), size(tt,2));
    figure; surf((CQ).^0.4,'edgecolor','none'); axis tight; view(0,90)
    title('Nearest neighbor interpolation with griddatan()');
    axis tight; view(0,90);
else
    disp('Data not interpolated. Not displaying final spectrogram.');
    CQ = 0;
end
