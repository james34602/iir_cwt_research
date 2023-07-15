lst_oct = round(linspace(16, 68, 8));
lst_fftLen = [4096];
numOvps = 4;
lst_hop = repmat(lst_fftLen, [numOvps, 1]) ./ (2.^(2 : (1 + numOvps)))';
lst_HFSamplingLimit = linspace(0.4, 1.3, 4);
coeffCollection = cell(length(lst_oct), length(lst_fftLen), size(lst_hop, 1), length(lst_HFSamplingLimit));
fs = 48000;
order = 2;
reqSynthesisWnd = 1;
delete(gcp('nocreate'));
parpool(2);
for lv1 = 1 : 8
    oct = lst_oct(lv1);
    for lv2 = 1 : length(lst_fftLen)
        fftLen = lst_fftLen(lv2);
        for lv3 = 1 : size(lst_hop, 1)
            hop = lst_hop(lv3, lv2);
            parfor lv4 = 1 : length(lst_HFSamplingLimit)
                HFSamplingLimit = lst_HFSamplingLimit(lv4);
                disp("Octave: " + string(oct) + ", FFT len: " + string(fftLen) + ", hop: " + string(hop) + ", hfSl: " + string(HFSamplingLimit))
                tic
                [coeff, f_q] = ltv_precomute2(fftLen, hop, fs, oct, order, HFSamplingLimit, 1, 1, reqSynthesisWnd, 2e-7);
                toc
                coeffCollection{lv1, lv2, lv3, lv4}.fftLen = fftLen;
                coeffCollection{lv1, lv2, lv3, lv4}.hop = hop;
                coeffCollection{lv1, lv2, lv3, lv4}.HFSamplingLimit = HFSamplingLimit;
                coeffCollection{lv1, lv2, lv3, lv4}.oct = oct;
                coeffCollection{lv1, lv2, lv3, lv4}.coeff = coeff;
                coeffCollection{lv1, lv2, lv3, lv4}.f_q = f_q;
            end
            clear coeff
            disp('Saving variables');
            save('matlab.mat')
        end
    end
end