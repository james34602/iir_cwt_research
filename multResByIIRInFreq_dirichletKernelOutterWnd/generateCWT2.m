lst_oct = round(linspace(20, 70, 8));
fftLen = 6144;
numOvps = 5;
lst_hop = repmat([fftLen], [numOvps, 1]) ./ (2.^(2 : (1 + numOvps)))';
coeffCollection = cell(length(lst_oct), size(lst_hop, 1));
fs = 48000;
order = 2;
reqSynthesisWnd = 1;
for lv1 = 1 : 8
    oct = lst_oct(lv1);
    for lv2 = 1 : size(lst_hop, 1)
        hop = lst_hop(lv2);
        HFSamplingLimit = 0.39;
        disp("Octave: " + string(oct) + ", FFT len: " + string(fftLen) + ", hop: " + string(hop) + ", hfSl: " + string(HFSamplingLimit))
        tic
        [coeff, f_q] = ltv_precomute2(fftLen, hop, fs, oct, order, HFSamplingLimit, 1, 1, reqSynthesisWnd, 5e-8);
        toc
        coeffCollection{lv1, lv2}.fftLen = fftLen;
        coeffCollection{lv1, lv2}.hop = hop;
        coeffCollection{lv1, lv2}.HFSamplingLimit = HFSamplingLimit;
        coeffCollection{lv1, lv2}.oct = oct;
        coeffCollection{lv1, lv2}.coeff = coeff;
        coeffCollection{lv1, lv2}.f_q = f_q;
    end
    clear coeff
    disp('Saving variables');
    save('matlab.mat')
end