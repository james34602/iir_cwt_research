# Frame based CWT

#### File List

1. multResByIIRInFreq_archived
   Initial prototype of the algorithm

2. multResByIIRInFreq_dirichletKernelOutterWnd
   The efficient frame-based CWT among all current projects

3. multResByIIRInFreq_generalizeOutterWnd
    Variant of the frame-based CWT that allows for a custom global analysis window

4. multResByIIRInFreq_FrameTheoretic
   Implementation of frame theoretic transform, provides a least-squares solution (canonical dual) for a nearly lossless inverse transform of the frame-based CWT

5. slidingCQT
   A CWT based on the Sliding DFT. Its main feature is outputting a spectral slice for each incoming sample, similar to a traditional CWT

6. multResByIIRInTime
   Implementation of aliased signal processing based CWT, it calculates the spectrum using a narrow Hann window, performs a circular shift on the frame, and applies IIR filtering to each frequency bin along the time axis to achieve multi-resolution analysis. This is an aliased signal processing method, and the input signal cannot(and not supposed) be reconstructed

7. minFunc_2012
   An unconstrained gradient-based optimizer

8. conv_constant-q-toolbox
   An inefficient, low-SNR, and invertible CQT implementation from previous work

9. discreteLevelWithWindowing
   Similar to multResByIIRInTime but without applying IIR filtering to each frequency bin along the time axis. The input signal cannot be reconstructed

10. iirLayer
    Derivative and test bench for time-varying State Variable Filter (SVF) in the time domain.

11. lsicqs
    An inefficient and complicated CQT implementation from previous work, targeted for specific applications

12. ADNode.m
    Reverse mode AD implementation in Matlab. Due to the fact that MATLAB's Symbolic Toolbox is too slow, I wrote an automatic (symbolic) differentiation calculator. For complicated formulas, ADNode is thousands or even tens of thousands of times faster than MATLAB's Symbolic Toolbox. It can solve for the reverse-mode AD derivatives of various functions, including the Discrete Fourier Transform, special cases of overlap-add, certain special functions, SVF filters (Backpropagation Through Time), smoothed round/floor/ceil functions, and anti-aliasing phase difference functions.

13. testMorletCWT.m
    A textbook implementation of the CWT. It includes two signal reconstruction (inverse transform) methods to demonstrate that as long as the signal is fully redundant, any filterbank can be perfectly reconstructed, including those with completely random filterbank coefficients

##### Pre-computed correction matrix

Link：[https://pan.baidu.com/s/14qUPAxsxLeBW1LjoNAea4A?pwd=uppl]()

##### Notes

1. I often refer to the CQT (Constant-Q Transform) as CWT. This should be acceptable, as all continuous-scale, multi-resolution transforms are forms of the CWT. The CQT is similar to the Morlet CWT, and the Q-factor parameter in CQT is analogous to the tunable time-bandwidth product parameter in the Morlet wavelet.
2. It took a month to write the derivative formulas and test code for the time-varying SVF filtering...

**Important Notice:** 

* **Commercial use of this project is strictly prohibited.**
*  **Redistribution of modified versions of this project without explicit permission is prohibited.**
