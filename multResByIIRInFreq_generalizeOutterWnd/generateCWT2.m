lst_oct = round(linspace(20, 70, 8));
lst_fftLen = [8192];
numOvps = 4;
lst_hop = repmat(lst_fftLen, [numOvps, 1]) ./ (2.^(2 : (1 + numOvps)))';
coeffCollection = cell(length(lst_oct), length(lst_fftLen), size(lst_hop, 1));
fs = 48000;
order = 2;
delete(gcp('nocreate'));
parpool(2);
for lv1 = 1 : 8
    oct = lst_oct(lv1);
    for lv2 = 1 : length(lst_fftLen)
        fftLen = lst_fftLen(lv2);
        parfor lv3 = 1 : size(lst_hop, 1)
            hop = lst_hop(lv3, lv2);
            HFSamplingLimit = 0.39;
            disp("Octave: " + string(oct) + ", FFT len: " + string(fftLen) + ", hop: " + string(hop) + ", hfSl: " + string(HFSamplingLimit))
            tic
            [coeff, f_q] = ltv_precomute2(fftLen, hop, fs, oct, order, HFSamplingLimit, 1, 1, 1e-7);
            toc
            coeffCollection{lv1, lv2, lv3}.fftLen = fftLen;
            coeffCollection{lv1, lv2, lv3}.hop = hop;
            coeffCollection{lv1, lv2, lv3}.HFSamplingLimit = HFSamplingLimit;
            coeffCollection{lv1, lv2, lv3}.oct = oct;
            coeffCollection{lv1, lv2, lv3}.coeff = coeff;
            coeffCollection{lv1, lv2, lv3}.f_q = f_q;
        end
        clear coeff
        disp('Saving variables');
        save('matlab.mat')
    end
end