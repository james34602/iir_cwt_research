%% STFChT
[spec, t, f, f0gram, f0s, f0_hyps_indxs, val_f0_hyps, selected_warps] = stfcht('pop1_long.wav');

%% PLOTS

A4 = 440; % reference A4 frequency for pitch grid in plots
chirps_from_labels = 0; % whether to plot chirp rates estimated from f0 labels

num_semitones_down = floor(12*log2(A4/f0s(1)));
num_semitones_up   = floor(12*log2(f0s(end)/A4));
value_ticks = A4*2.^((-num_semitones_down:num_semitones_up)/12);
pos_ticks = interp1(f0s, 1:length(f0s), value_ticks);
f_ylabels = num2str(value_ticks', '%4.2f');

notes_names = {'A-'; 'Bb'; 'B-'; 'C-'; 'C#'; 'D-'; 'Eb'; 'E-';  'F-'; 'F#'; 'G-'; 'Ab'} ;
semitones = (-num_semitones_down:num_semitones_up)';
octaves = num2str(floor((semitones+9)/12) + 4,'%1d');
semitones_notes_names = notes_names(mod(semitones,12)+1);
ylabels = semitones_notes_names;
for i=1:length(semitones)
    ylabels{i} = strcat(ylabels{i}, octaves(i));
    ylabels{i} = strcat(ylabels{i}, '-');
    ylabels{i} = strcat(ylabels{i}, f_ylabels(i,:));
end

% plot of Short Time Fan Chirp Transform
figure(1)
spec = spec ./ max(abs(spec(:)));
imagesc(t,f,20*log10(abs(spec)))
caxis([-90, 0])
set(gca,'YDir','normal');
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Short Time Fan Chirp Transform')

% plot of F0gram
figure(2)
hh = imagesc(t,1:length(f0s),double(f0gram));
set(gca,'YDir','normal');
title('F0gram for the best direction for each f0')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(get(hh,'Parent'),'YTick',pos_ticks)
set(get(hh,'Parent'),'YTickLabel',ylabels)
ylim([1 length(f0s)])
grid

% same plot with f0 candidates
hop = 3;
f0_hyps = f0_hyps_indxs(:,1:hop:end);
t_f0_hyps = t(1:hop:end);

figure(3)
hh = imagesc(t,1:length(f0s),double(f0gram));
set(gca,'YDir','normal');
title('F0gram and f0 candidates')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(get(hh,'Parent'),'YTick',pos_ticks)
set(get(hh,'Parent'),'YTickLabel',ylabels)
hold on
plot(t_f0_hyps,f0_hyps(1,:),'ok','MarkerSize',4)
plot(t_f0_hyps,f0_hyps(2,:),'+k','MarkerSize',4,'LineWidth',0.1)
plot(t_f0_hyps,f0_hyps(3,:),'xk','MarkerSize',4,'LineWidth',0.1)
legend('First','Second','Third');
hold off
ylim([1 length(f0s)])
grid


% selected chirp rates
figure(4)
plot(t,selected_warps(1,:),'ok', 'MarkerSize',4);
if chirps_from_labels == 1
    hold on;
    plot(tchr,chirp_rates,'-xk','Color',[0.6,0.6,0.6], 'MarkerSize',4);
    hold off
    legend('selected chirp rates', 'chirp rates from labels')
end
xlabel('Time (s)')
ylabel('\alpha')
title('Chirp rate estimation for best candidate')
axis tight
grid
function [f0gram t ind_maxs f0_hyps_indxs val_f0_hyps spec selected_warps f0_medians chirp_rates tchr] = compute_fcht(y, fs, nfft, fmax, hop, warps, f0s, accums, f0_params, glogs_params)
% computes a spectrogram based on the FChT, a f0gram and predominant f0 hypotesis
%
%  in:
%               y - audio signal
%            nfft - number of fft points
%            fmax - maximum frequency of interest (Hz)
%             hop - hop in samples
%
%           warps - warpings parameters design
%                   data structure to do the interpolation using "c" code
%                   with fields:
%                    .pos1.int       - "c" index of previous sample
%                    .pos1.frac      - fractional value
%                    .fs_orig        - sampling frequency after oversampling
%                    .fs_warp        - sampling frequency of warped signal
%                    .nsamps_torig   - number of samples in the oversampled signal frame
%                    .fact_over_samp - oversampling factor
%                    .chirp_rates    - chirp rate values
%
%             f0s - f0 grid values where to compute the f0gram
%
%          accums - spectrum array structs for efficient computation of the
%                   log spectrum accumulation, efficient attenuation of
%                   mutliples and submultiples, and high pass filtering
%
%                   .harmonic_accumulations: an array of structs with fields,
%                           .pos_int:       array of integer indexes for fast interpolation
%                           .pos_frac:      array of fractionary part for fast interpolation
%                           .weights:       array of weights to be aplied to addends
%                           .accum_length:  array of lengths of accumulations
%                           .global_weight: array of global weights applied to the accumulations
%
%                   .harmonic_attenuations: struct with the same fields as above
%
%                   .HP_logS: if high pass filtering is applied (1:yes 0:no)
%                   .HP_logS_components: values used in high pass filtering
%
%       f0_params - f0-gram parameters, struct with fields:
%                    .f0_params.f0min           - minimun fundamental frequency
%                    .f0_params.num_octs        - number of octaves of the f0-gram
%                    .f0_params.num_f0s_per_oct - number of f0s per octave
%                    .f0_params.num_f0_hyps     - number of f0 hypotesis to extract
%                    .f0gram_type               - type of f0gram to return:
%                                                  'none' - no f0gram is returned (default)
%                                                'single' - single optimum direction f0gram
%                                              'compound' - best direction for each f0
%
%    glogs_params - gathered log spectrum parameters, struct with fields:
%                    .HP_logS             - whether to high-pass the GlogS
%                    .att_subharms        - whether to attenuate subharmonics
%                    .apply_GlogS_model   - whether to apply the normalization model
%                    .compute_GlogS_model - whether to learn the normalizartion model
%                    .median_poly_coefs   - model parameter value (median correction)
%                    .sigma_poly_coefs    - model parameter value (sigma  correction)
%
%   labels_params - parameters related to use of labels of fundamental
%                   frequency as the ones provided in the RWC or MIREX
%                   melody extraction test set, struc with fields:
%                   .chirp_rate_from_labels - use best chirp rate according to f0 labels
%                   .only_melody_frames     - process only frames labeled as melody
%                   .database               - name of the database ('mirex','rwc')
%                   .database_path          - path to the database files
%                   .audio_file             - audio filename (used to set labels filename)
%                   .fs_database            - sampling frequency of the database files
%                   .path_mat               - path to save mat files with the computations
%
% out:
%
%          f0gram - F0gram, salience value for each f0 of the grid along audio frame instants
%               t - central time instants of frames
%        ind_maxs - warping index of optimal warping for each frame
%   f0_hyps_indxs - indexes of the f0 candidates in the grid for each frame
%     val_f0_hyps - F0gram amplitude value for each f0 candidate for each frame
%            spec - spectrum for optimal warping for each frame, i.e. STFChT
%  selected_warps - indexes of selected warps for each f0 hypotesis for each frame
%
% if chirp rates are estimated from f0 labels (labels_params.chirp_rate_from_labels > 0):
%      f0_medians - median f0 estimated for each frame using labels
%     chirp_rates - chirp_rate estimated from labels for each frame
%            tchr - frame time instants of estimations from labels
%
% CWT
oct = 90;
reqSynthesisWnd = 1;
% number of samples of the warped signal frame
fftLen = size(warps.pos1_int,1);
halfLen = getFFTHalfLen(fftLen);
% number of points of pre and post padding used to set initial conditions
thetas = 0:(fftLen/2);
thetas(1) = eps;
%% Compute poles
sigmas = thetas / (oct * pi);
[b, a, c1, c2] = gauss_precompute(sigmas);
%% Pole limiting
analysisWnd = hann(fftLen, 'periodic') .^ (1 / reqSynthesisWnd);
synWnd1 = hann(fftLen, 'periodic') .^ reqSynthesisWnd;
chopedWnd1 = analysisWnd(halfLen : end);
chopedWnd2 = synWnd1(halfLen : end);
halfWndLen = halfLen - 1;
digw2 = linspace(0, pi, halfLen);
digw = digw2(1 : halfWndLen);
cplxFreq = exp(1i*digw); % Digital frequency must be used for this calculation
h = (cplxFreq .* cplxFreq .* b(:)) ./ (cplxFreq .* (cplxFreq + a(:, 2)) + a(:, 3));
h2 = (h .* conj(h)) .* chopedWnd1.';
h2 = h2 .* chopedWnd2.';
theoreticalWindowShape = [zeros(size(thetas, 2), 1), h2(:, (halfLen-1):-1:2), h2];
tol = min(1.5, max(1, (fftLen / hop) / 5));
hopsizeTol = min(fftLen, ceil(hop * 2 * tol));
if mod(hopsizeTol, 2) == 1
    hopsizeTol = hopsizeTol - 1;
end
smallestPossibleWnd = [zeros((fftLen - hopsizeTol) / 2, 1); hann(hopsizeTol, 'periodic'); zeros((fftLen - hopsizeTol) / 2, 1)];
wndDif = theoreticalWindowShape' - smallestPossibleWnd;
wndDifPwr = sum(abs(wndDif), 1);
[~, firstUndersampling1] = min(wndDifPwr);
firstUndersampling2 = find(any(wndDif < 0), 1, 'first');
firstUndersampling = ceil((firstUndersampling1 + firstUndersampling2) / 2);
% firstUndersampling = max(firstUndersampling, fix(fftLen / hop * oct / 2));
if ~isempty(firstUndersampling)
    thetas = 0:(fftLen/2);
    thetas(1) = eps;
    thetas(firstUndersampling : end) = thetas(firstUndersampling);
    thetas1 = thetas;
    thetas1(halfLen+1:fftLen) = conj(thetas1(halfLen-1:-1:2));
else
    thetas1 = thetas;
    thetas1(halfLen+1:fftLen) = conj(thetas1(halfLen-1:-1:2));
end
% Eliminate oscillation around corner
time = 0.026 * fftLen; % More elements in array less smoothing is needed
alpha = 1 / (1 + time);
[b2, a2] = butter(1, alpha);
thetas1 = filtfilt(b2, a2, thetas1);
mtheta = 10;
if min(thetas1) > mtheta
    flr = min(thetas1);
    cel = max(thetas1);
    thetas1 = thetas1 - flr;
    thetas1 = thetas1 ./ max(thetas1);
    thetas1 = thetas1 * (cel - mtheta) + mtheta;
end
thetas1(thetas1 < eps) = eps;
sigmas = (thetas1 ./ fftLen) ./ oct / pi * fftLen;
[b, a, c1, c2] = gauss_precompute(sigmas);
% number of warpings
num_warps = size(warps.pos1_int,2);
% if chirp rate is derived from labels

% windows applied to warped signal frames
windows = repmat(analysisWnd,[1,num_warps]);

% frame indexes
frame_indexes = (1:warps.nsamps_torig);
% time sample indexes of each signal frame
frame_time_samples = (1:(warps.fact_over_samp * hop):warps.fact_over_samp*(length(y)-warps.nsamps_torig));
% number of frames
nframes = length(frame_time_samples);
% number of f0 values
nf0s = length(f0s);
% frequency bin index corresponding to fmax
ind_fmax = ceil(fmax/warps.fs_warp*nfft);
% minimun relative distance between f0 hypotesis peaks
min_peaks_dist = 1/100;

% optimum f0gram
if ~strcmp(f0_params.f0gram_type, 'none')
    f0gram = zeros(nf0s, nframes, 'single');
else
    f0gram = [];
end
% warping index of optimal warping
ind_maxs = zeros(nframes,1);
% indexes of f0 hypotesis for each frame
f0_hyps_indxs = zeros(f0_params.num_f0_hyps,nframes);
% values of f0 hypotesis for each frame
val_f0_hyps = zeros(f0_params.num_f0_hyps,nframes);
% indexes of selected warps for each f0 hypotesis for each frame
selected_warps = zeros(f0_params.num_f0_hyps,nframes);
% spectrum for optimal warping
spec = zeros((ind_fmax+1), nframes, 'single');

% ====== GLOGS CORRECTION MODEL ======
if glogs_params.apply_GlogS_model
    % fundamental frequency indexes
    pos = 1:length(f0s);
    % model correction
    glogs_params.median_correction = polyval(glogs_params.median_poly_coefs,pos)';
    glogs_params.sigma_correction  = polyval(glogs_params.sigma_poly_coefs, pos)'*polyval(glogs_params.sigma_poly_coefs,pos(1));
end

% ====== PREPROCESSING ======
% oversampling
y = resample(y,warps.fact_over_samp,1);

% ======  PROCESSING  ======
% frame time instants
t = (frame_time_samples + warps.nsamps_torig/2)/warps.fs_orig;

index_of_frames = 1:1:nframes;

% compute f0 preference function, if needed
if f0_params.prefer == 1
    f0_preference_weights = repmat(f0_preference_weighting(f0s, f0_params),[1,num_warps]);
end

for pos = index_of_frames
    % signal frame
    i = frame_time_samples(pos);
    x_frame = y(i+frame_indexes,1);

    index_warps = 1:num_warps;

    % warped frames
    x_warpings = quick_interp1_int_frac_c(x_frame, warps.pos1_int(:,index_warps), warps.pos1_frac(:,index_warps));

    % fft of warped frames (FChT)
    fcht = fft(fftshift(x_warpings .* windows, 1));
    tmp = zeros(size(fcht, 1), 1, 'like', fcht);
    for alps = 1 : size(x_warpings, 2)
        fcht(:, alps) = ltv1Slice(fcht(:, alps), tmp, b, a, c1, c2) ./ 1;
    end

    % consider spectrums only up to fmax
    fcht = fcht(1:(ind_fmax+1),:);

    % log spectrum
    logS = log10_1_alpha_abs_c(fcht,10); % logS = log10(1+10*abs(fcht));

    % harmonic accumulation to compute pitch salience
    f0gram_frame = harmonic_accumulation(num_warps, logS, accums.harmonic_accumulations, accums.harmonic_attenuations, nf0s, f0_params.num_f0s_per_oct, accums.HP_logS_components, glogs_params);

    % apply f0 preference weighting, if needed
    if f0_params.prefer == 1
        f0gram_frame = f0gram_frame .* f0_preference_weights;
    end

    % predominant pitch estimation
    [~, indx] = max(max(f0gram_frame));
    ind_maxs(pos) = indx;

    % save optimum spectrum (if needed) and f0gram
    spec(:,pos) = fcht(:,indx);

    % whether to return a single optimum direction f0gram
    if strcmp(f0_params.f0gram_type, 'single')
        f0gram(:,pos) = f0gram_frame(:,indx);
    end
    % whether to compute compound f0gram and return it
    if strcmp(f0_params.f0gram_type, 'compound')
        % optimum warping for each f0
        [f0gram_maxs, chirp_rate_indexes] = max(f0gram_frame,[],2);
        if strcmp(f0_params.f0gram_type, 'compound')
            f0gram(:,pos) = f0gram_maxs;
        end
    end

    % compute multiple f0 hypotesis
    % find local maxima of compound f0gram
    ind_maxs = [];
    for j = 2:nf0s-1
        if ((f0gram_maxs(j) > f0gram_maxs(j-1)) && (f0gram_maxs(j) > f0gram_maxs(j+1)))
            ind_maxs = [ind_maxs; j];
        end
    end
    % sort local maxima
    [~, inds] = sort(f0gram_maxs(ind_maxs),'descend');
    % save biggest local maxima
    num_max_inds = length(inds);
    if num_max_inds > 0
        num_inds = min(f0_params.num_f0_hyps,num_max_inds);
        % add f0 candidates considering they are 1% appart
        f0_hyps_indxs(1,pos) = ind_maxs(inds(1));
        val_f0_hyps(1,pos) = f0gram_maxs(f0_hyps_indxs(1,pos));
        selected_warps(1,pos) = warps.chirp_rates(chirp_rate_indexes(ind_maxs(inds(1))));
        k = 2;
        for m = 2:num_inds
            fin = 0;
            while ((k <= num_max_inds) && ~fin)
                fin = (min(abs(f0s(f0_hyps_indxs(1:m-1,pos)) - f0s(ind_maxs(inds(k))))./f0s(f0_hyps_indxs(1:m-1,pos))) > min_peaks_dist);
                if ~fin
                    k = k + 1;
                end
            end
            if fin
                f0_hyps_indxs(m,pos) = ind_maxs(inds(k));
                val_f0_hyps(m,pos) = f0gram_maxs(f0_hyps_indxs(m,pos));
                selected_warps(m,pos) = warps.chirp_rates(chirp_rate_indexes(ind_maxs(inds(k))));
                k = k + 1;
            end
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function accum = accumulate_vectors(x,accum_interpolations,num_f0s_per_oct)
% function accum = accumulate_vectors(x,accum_interpolations)
%
% accumulation of the signal x in position given in accum_interpolations
%
% this is done to compute the salience of diferent f0 candidates, by
% accumulation of the spectrum at the position of harmonics
% each accumulation is then weighted by a global weigth that in this case
% corresponds to the number of accumulated points
%

% interpolation of values at postions given by accum_interpolations
interps = quick_interp1_int_frac_c(x,accum_interpolations.pos_int,accum_interpolations.pos_frac);

%accum = accum_arrays(interps,accum_interpolations.accum_length);
accum = recursive_octave_accum(accum_arrays(interps,accum_interpolations.accum_length)',num_f0s_per_oct)';

% weighting of accumulated values by the number of accumulated points
accum = accum.*accum_interpolations.global_weight;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f0gram_frame = harmonic_accumulation(num_warps, logS, accum_interpolations, subharmonic_interpolations, num_f0s, num_f0s_per_oct, HP_logS_components, glogs_params)
%function f0gram_frame = harmonic_accumulation(num_warps, logS, accum_interpolations, subharmonic_interpolations, num_f0s, num_f0s_per_oct, HP_logS_components, HP_logS, att_subharms, compute_GlogS_model, apply_GlogS_model, median_correction, sigma_correction)

% high pass filter and mean substraction
if(~isempty(HP_logS_components))
    if(glogs_params.HP_logS == 1)
        logS = logS - HP_logS_components*(logS'*HP_logS_components)';
    end
    for j = 1:num_warps
        if(glogs_params.HP_logS == 2)
            logS(:,j) = quick_filtfilt_1_v2010(B,A,logS(:,j));
        end
        logS(:,j) = (logS(:,j)-mean(logS(:,j)));
    end
end

% global variables to learn the GlogS model
if glogs_params.compute_GlogS_model
    global mean_glogs
    global mean_glogs_post
    global mean_glogs2
    global mean_glogs_post2
    global cant_N
end

f0gram_frame = zeros(num_f0s, num_warps);

for j = 1:num_warps

    % accumulation of log spectrum
    glogs = accumulate_vectors(logS(:,j),accum_interpolations, num_f0s_per_oct);

    % post-processing of glogs for subharmonic attenuation
    f0gram_frame(:,j) = post_process_glogs(glogs, subharmonic_interpolations, num_f0s, glogs_params.att_subharms, num_f0s_per_oct);

    % apply model correction
    if glogs_params.apply_GlogS_model
        %f0gram_frame(:,j) = (f0gram_frame(:,j)-polyval(P,pos)')./(polyval(Ps,pos)'*polyval(Ps,pos(1)));
        f0gram_frame(:,j) = (f0gram_frame(:,j)-glogs_params.median_correction)./glogs_params.sigma_correction;
    end

    % model learning
    if glogs_params.compute_GlogS_model
        if(isempty(mean_glogs))
            mean_glogs = glogs;
            mean_glogs_post = f0gram_frame(:,j);
            mean_glogs2 = glogs.*glogs;
            mean_glogs_post2 = f0gram_frame(:,j).*f0gram_frame(:,j);
            cant_N = 1;
        else
            cant_N = cant_N + 1;
            mean_glogs = (mean_glogs*(cant_N-1) + glogs)/cant_N;
            mean_glogs_post = (mean_glogs_post*(cant_N-1) + f0gram_frame(:,j))/cant_N;
            mean_glogs2 = (mean_glogs2*(cant_N-1) + glogs.*glogs)/cant_N;
            mean_glogs_post2 = (mean_glogs_post2*(cant_N-1) + f0gram_frame(:,j).*f0gram_frame(:,j))/cant_N;
        end
    end

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glogs_pospro = post_process_glogs(glogs, subharmonic_interpolations, num_f0s, att_subharms, num_f0s_per_oct)
% function glogs_pospro = post_process_glogs(glogs, subharmonic_interpolations, num_f0s, att_subharms, num_f0s_per_oct)
%
% postprocessing of the gathered log spectrum by attenuation of each value
% by the maximum of the values at the subharmonics
% and subharmonics attenuations by the scaled value at the next octave

if ~att_subharms
    glogs_pospro = glogs(end-num_f0s+1:end) - max_arrays(quick_interp1_int_frac_c(glogs,subharmonic_interpolations.pos_int,subharmonic_interpolations.pos_frac), subharmonic_interpolations.accum_length);
else
    glogs_pospro = glogs((num_f0s_per_oct+1):end) - max_arrays(quick_interp1_int_frac_c(glogs,subharmonic_interpolations.pos_int,subharmonic_interpolations.pos_frac), subharmonic_interpolations.accum_length);
    %glogs_pospro = glogs_pospro(1:num_f0s) - 1/2 * glogs_pospro(num_f0s_per_oct+1:end);
    glogs_pospro = glogs_pospro(1:num_f0s) - 1/3 * glogs_pospro(num_f0s_per_oct+1:end);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f0_preference_weights = f0_preference_weighting(f0s, f0_params)
% models f0 preference, for instance for melody detection, as a gaussian
% whose parameters are derived for example from labeled data

MIDI_grid = 69 + 12 * log2(f0s/440);

bell = 1/sqrt(2*pi*f0_params.prefer_stdev^2).*exp(-(MIDI_grid-f0_params.prefer_mean).^2/(2*f0_params.prefer_stdev^2));

% kind of simple smoothing for non observed frequencies
smoothing_offset = 0.01;
f0_preference_weights = ( smoothing_offset + bell') / ( smoothing_offset + 1);

end

function [spec, t, f, f0gram, f0s, f0_hyps_indxs, val_f0_hyps, selected_warps] = stfcht(audio_file)
% read audio file
[y, fs] = audioread(audio_file);
% script to set FChT analysis parameters
% see below for a description of each parameter
%% =============  WARPING PARAMETERS  =============
% maximum frequency of interest
fmax = fs / 2; % Hz
% number of samples of the warped signal frame
warp_params.nsamps_twarp = 4096;
% number of fft points (controls zero-padding)
% maximum value of normalized frequency deviation (alpha)
warp_params.alpha_max = 6;
% number of warpings
warp_params.num_warps = 21;
% oversampling factor
warp_params.fact_over_samp = 4;
% hop in samples
hop = 256;

%% =============   F0-GRAM PARAMETERS  =============
% minimun fundamental frequency
f0_params.f0min = 80; % Hz
% number of octaves
f0_params.num_octs = 4;
% number of f0s per octave
f0_params.num_f0s_per_oct = 192;
% number of f0 hypotesis to extract
f0_params.num_f0_hyps = 10;
% type of f0gram returned: 'none', 'single', 'compound'
f0_params.f0gram_type = 'compound';
% whether to use a f0 preference guassian function
f0_params.prefer = 1;
% mean of f0 preference guassian function (MIDI number for C4)
f0_params.prefer_mean = 60;
% stdev of f0 preference guassian function (stdev in MIDI numbers)
f0_params.prefer_stdev = 18;

%% ======== GATHERED LOG SPECTRUM PARAMETERS =======
% high-pass logS
glogs_params.HP_logS = 1;
% whether to attenuate subharmonics
glogs_params.att_subharms = 1;
% whether to apply the GlogS correction model
glogs_params.apply_GlogS_model = 1; % may be changed later
% whether to learn the GlogS correction model from a labeled database
glogs_params.compute_GlogS_model = 0; % may be changed later
% model parameter variables (default values)
glogs_params.median_poly_coefs = [-0.000000058551680 -0.000006945207775 0.002357223226588];
glogs_params.sigma_poly_coefs  = [ 0.000000092782308  0.000057283574898 0.022199903714288];

% add path to the c_code and m_code directories
% addpath(['.' dir_slash 'c_code']);
% addpath(['.' dir_slash 'm_code']);

% design for efficient fcht computation
[warps f0s accums] = design_fcht(fs, warp_params.nsamps_twarp, fmax, warp_params, f0_params, glogs_params);

% fcht computation
tic
[f0gram t ind_maxs f0_hyps_indxs val_f0_hyps spec selected_warps] = compute_fcht(y, fs, warp_params.nsamps_twarp, fmax, hop, warps, f0s, accums, f0_params, glogs_params);
toc

ind_fmax = ceil(fmax / warps.fs_warp * warp_params.nsamps_twarp);
f = linspace(0,fmax,ind_fmax+1);
end
function halfLen = getFFTHalfLen(fftLen)
if mod(fftLen, 2) == 0
    halfLen = (fftLen / 2) + 1;
else
    halfLen = (fftLen + 1) / 2;
end
end