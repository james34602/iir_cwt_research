What is it?
    -A fully automated flowgraph (signal flow graph) analysis tool using 
     Matlab’s symbolic and control system toolbox. 
    -A tool for generating one or several transfer functions for a given 
     control system.
    -An analysis tool for both continuous- and discrete-time control systems. 
    -Includes a function for finding a system of equations

Introduction
    From a simple textual nodelist model the user can easily generate one or 
    several transfer functions between any nodes and apply the result 
    to control system analysis. Such a versatile tool can be applied to 
    filter design, delta-sigma modulator analysis or even reflection 
    analysis s-parameter measurements in RF techniques and in optics. 
    Some design examples are attached to the package and also presented in the manual.

Instructions
    -Unpack ZIP-file and extract the file examples.tar in the same directory.
    -Read the manual.
    -Study the demo (flowdemo.m)
    -See what's in the *.flw files.
    -See help for flow_tf.m.
    -Output transfer function(s) are struct arrays,
     so typing e.g. H.sym{1} shows the result (see the demo).

Contents
    -m-files (*.m):
        flow_tf:            the flowgraph analysis function
        find_soe:           find the system of equations
        flowdemo:           a demo script
        test_all_nodelists: runs all *.flw files in the package

    -textual nodelist files (in examples.tar):
        zpar            a terminated z-parameter two-port circuit model
        biquad:         biquad filter example
        IIR:            IIR filter example
        integra:        a simple integrator
        efb:            error-feedback delta sigma
        MOD1:           1st order delta sigma
        MODL:           Lth order pure differentiation delta sigma
        CRFBSYM:        3rd order symbolic delta-sigma
        CRFBNUM:        3rd order numeric delta-sigma
	    fiori           3rd order efb delta-sigma. 
                        From paper by Fiori&Maloberti at ecctd2005

    -manual.pdf

    -readmefirst.txt

Authors:
  Marko Neitola, Timo Rahkonen 
  Electronics laboratory
  Dept. of Electrical and Information Engineering
  University of Oulu, Finland. 