while 1
    ckmin = randi(10, 1, 1);
    fftLen = randi(100, 1, 1);
    halfLen = randi(100, 1, 1);
    if halfLen > fftLen
        halfLen = fftLen;
    end
    halfLen = fftLen;
    ck = max(ckmin,floor(fftLen/halfLen)); % determines size of matrix Phi (and reconstruction error)
    Nx = ck*fftLen;
    Nm = ceil(Nx/hop);
    k = ck*ceil(fftLen/hop);
    if Nm ~= k
        disp([ckmin, fftLen, halfLen]);
    end
end