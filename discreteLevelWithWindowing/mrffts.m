o = mrfft(data,fs);
magSpec = abs(o)';
imagesc(power(magSpec, 0.1))
colormap(jet);
set(gca,'YDir','normal');
function Pmr = mrfft(x,Fs)
% MRFFT   Produce and display a multi-resolution spectrogram
%   MRFFT(X,N,L,FS,FMAX) takes waveform data X sampled at FS and calculates the N-point
%   FFT at varying resolutions specified by buckets of the Bark scale up to FMAX given in Hz.
%   The number of different resolutions is determined by L, the hop size,
%   where N/L = 8.
%
%   For now, N and L are fixed at 2048 and 256, respectively
%   fmax must also be > 3150 Hz
%
%   Bark scale:
%   0-630 Hz: highest freq res (i.e. # windows/N points = r = 8)
%   630-1480 Hz: r = 4
%   1480-3150 Hz: r = 2
%   3150-fmax: r = 1
%

h = [0.53836, 0.46164]; % Hamming
h = [0.5 0.5]; % Hann
N = 4096;
L = 512;
rat = log(N / L) / log(2);

halfLen = N / 2 + 1;
nframes = floor(length(x)/L);
% circular buffer for storing elementary FFTs
% need additional N/L slots due to shifts during freq windowing
buffer = zeros(N/L,halfLen+N/L);
Pmr = zeros(nframes, halfLen);
for l = 0:nframes-1
    c = mod(l,N/L);
    xc = zeros(N,1);
    xc(c*L+1:(c+1)*L) = x(1+l*L:(l+1)*L); % Zero padded hoping segment
    Xc = fft(xc,N); % DFT of zero padded hoping segment
    buffer(c+1,:) = Xc(1:halfLen+N/L);
    for r = 2.^(0:rat)
        M = r*L;
        Xr = zeros(1,halfLen+N/L);
        cindex = c+1;
        if (l >= r-1)
            for i = 1:r
                Xr = Xr + buffer(cindex,:);
                if (cindex == 1)
                    cindex = N/L;
                else
                    cindex = cindex - 1;
                end
            end
            Xr = Xr .* exp(1i*2*pi*(0:halfLen+N/L-1)*(c-r+1)*L/N); % "twiddle" factor

            % frequency windowing
            Xr = h(1)*Xr(1:halfLen) - h(2)/2*([conj(fliplr(Xr(2:r*N/M+1))) Xr(1:halfLen-r*N/M)] + Xr(r*N/M+1:end));
            if (r==8)
                krange = 1:floor(630/Fs*N);
            elseif (r==4)
                krange = floor(630/Fs*N)+1:floor(1480/Fs*N);
            elseif (r==2)
                krange = floor(1480/Fs*N)+1:floor(3150/Fs*N);
            else
                krange = floor(3150/Fs*N)+1:halfLen;
            end
            Pmr(l+1, krange) = abs(Xr(krange))';
        end
    end
end
end