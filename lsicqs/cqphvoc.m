% %%%%%%%%%%%%%%%%%%%%%%%%% Version Dec 17, 2011 %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Uses the phase vocoder resysnthesis idea to estimate 'fake' phase values
% for the inverse magnitude transform of the constant Q transform.
% The CQ transform is assumed to have been produced by using the transform
% kernel matrix generated using cqker4.m, with the phase values discarded.
% cf. Phase vocoder PV.m Bill Sethares 5/2005
%
% CQold = Given CQ spectrogram
% CQmagph = CQ spectrogram with phase values reconstructed
% timeMap = corresponding time values for each element in the CQ specgram
% f = vector of frequencies that were analyzed in the CQ spgram
% sr = sampling rate
% bpo = bins per octave
%
% Keep the "time streched" timeMap ready
% for reference and use those time values directly to compute phase deltas.
%
% Mar 15 2011: Bug fix -- assignedtime changed to assignedtime_outslices in
% line 189.
%
% Dec 17, 2011: Added documentation comments throughout the code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CQmagph = cqphvoc(z, f, sr, bpo, op, timeMap, timeMap1p5, Nvec, kk, scale)
% z is the sound signal here, CQ is the vectorized complex CQSTFT
% CQold is the given complex spectrogram (reshaped from CQ)

%freqidxnumber = 30;
analysistimeslicenumber = length(f);

%t=0:1/sr:10; z = sin(2*pi*f(freqidxnumber)*t);
%z = z(1:1116)';

%z = wavread('srgm.wav');
%z = resample(z, 4, 147);
%z = z(1:1116);

%zstruct = load('gong');
%z = resample(zstruct.y, 75, 512);
%z = z(1:1116);

CQ = cell2mat(kk)*z;
CQtemp = CQ;
CQold = cell(size(kk));
for i=1:length(CQold)
    CQold{i} = (CQ(1:size(kk{i},1)));
    CQ = CQ(size(kk{i},1)+1:end);
end
CQ = CQtemp;
clear CQtemp


MAX_PEAK = 5;
EPS_PEAK = 0.5;
all2pi = 2*pi*(0:100);
Q = 1/(2^(1/bpo) - 1);
jp = 100-op;

CQmagph = (timeMap1p5); % we want both to have identical structure
for i =1:length(CQmagph)
    if length(CQmagph{i})<=length(CQold{i})
        CQmagph{i} = CQold{i}(1:length(CQmagph{i}));
    else
        CQmagph{i} = [CQold{i}; CQold{i}(end)*ones(length(CQmagph{i})-length(CQold{i}),1)];
    end
end

assignedtime_outslices = (CQmagph);
for i=1:length(assignedtime_outslices)
    assignedtime_outslices{i} = NaN.*ones(size(assignedtime_outslices{i}));
end % NaN indicates a particular point has never been assigned a phase

dtin = (Nvec(analysistimeslicenumber)*jp/100)/(2*sr);
times = timeMap{analysistimeslicenumber}(1):dtin:timeMap{analysistimeslicenumber}(end)+dtin;

% Set output hop to a multiple of input hop for time-scale modification
dtout = scale*dtin;

times_outslices = times(1):dtout:(dtout/dtin+1)*times(end);
times_outslices = times_outslices(1:length(times));

phadvance = zeros(size(timeMap));
ph = zeros(size(timeMap));
phold = ph;

for k=1:length(times)   % go to each timeslice
    curtime = times(k);
    mag = zeros(size(timeMap));
    idxvec = zeros(size(timeMap));
    tdeltavec = zeros(size(timeMap));
    for i=1:length(timeMap)
        tmp = ( timeMap{i} - curtime );
        abstmp = abs(tmp);
        [~, idx] = min(abstmp);
        mag(i) = abs(CQold{i}(idx));
        idxvec(i) = idx;
        tdeltavec(i) = tmp(idx);
    end
    
    curtime_outslices = times_outslices(k);
    idxvec_outslices = zeros(size(timeMap1p5));
    tdeltavec_outslices = zeros(size(timeMap1p5));
    for i=1:length(timeMap1p5)
        tmp = ( timeMap1p5{i} - curtime_outslices );
        abstmp = abs(tmp);
        [~, idx] = min(abstmp);
        idxvec_outslices(i) = idx;
        tdeltavec_outslices(i) = tmp(idx);
    end
    
    for i=1:length(ph)
        ph(i) = angle(CQold{i}(idxvec(i))) - 2.*pi.*f(i).*tdeltavec(i); 
    end
    
    % See the phase vocoder implementation by Sethares
    % (http://sethares.engr.wisc.edu/vocoders/phasevocoder.html)    
    peaks = findPeaks4(mag, MAX_PEAK, EPS_PEAK, 0); % gives magnitude peaks and regions, smallest index to biggest index
    [~, inds] = sort(mag(peaks(:,2))); % sorts magnitude peaks by height, small to big and gets the index ordering for that
    peaksort = peaks(inds,:); % reorder peaks in increasing order of height of the peaks
    pc = peaksort(:,2); % get peak locations in increasing order of size of peak
    
    % Frequency estimation step (using Eq. (6) in the LSICQS paper)
    bestf = zeros(size(pc));
    for tk=1:length(pc)
        dtheta = (ph(pc(tk)) - phold(pc(tk)))+all2pi;
        fest = dtheta/(2*pi*dtin);                     % analyze with dtin
        [~, indf] = min(abs(f(pc(tk))-fest));
        bestf(tk) = fest(indf);
    end
    
    % 
    magout = mag; phout = ph;
    for tk = 1:length(pc)
        fdes = bestf(tk);
        freqind = (peaksort(tk,1):peaksort(tk,3)); % Extract the peak region
        magout(freqind) = mag(freqind);
        % Apply time-slice correction for phase at the peak (Fig (1) in the
        % LSICQS paper
        phadvance(peaksort(tk,2)) = phadvance(peaksort(tk,2))+...
            2*pi*fdes*dtout;
        % Set phase in peak-region equal to that at peak by Theorem (1) in
        % the LSICQS paper
        phout(freqind) = ones(size(freqind)).*phadvance(peaksort(tk,2));        
    end
    % correct for timeslice phase when referenced at time slice in new
    % timeMap, again refer Fig (1) in LSICQS paper
    phcorrection = 2.*pi.*f'.*tdeltavec_outslices;
    phout = phout+phcorrection;
    phout = wrapToPi(phout);
    % assign these phases into CQmagph
    for i=1:length(timeMap1p5)
        [re, im] = pol2cart(phout(i), abs( CQmagph{i}(idxvec_outslices(i))) );
        if isnan(assignedtime_outslices{i}(idxvec_outslices(i)))
            assignedtime_outslices{i}(idxvec_outslices(i)) = curtime_outslices;
            CQmagph{i}(idxvec_outslices(i)) = re + 1i*im;
        elseif abs(timeMap1p5{i}(idxvec_outslices(i)) - assignedtime_outslices{i}(idxvec_outslices(i)))>...
                abs(timeMap1p5{i}(idxvec_outslices(i)) - curtime_outslices)
            assignedtime_outslices{i}(idxvec_outslices(i)) = curtime_outslices;
            CQmagph{i}(idxvec_outslices(i)) = re + 1i*im;
        end
    end
    phold = ph;
    
end