if ismac || isunix
    mex -R2018a -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -DNOCHECK -output ltv ltv.c
    mex -R2018a -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -DNOCHECK -output ltv1Slice ltv1Slice.c
    mex -R2018a -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -DNOCHECK -output cplxReassignment cplxReassignment.c
    mex -R2018a -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -DNOCHECK -output fcwtComplexity fcwtComplexity.c
elseif ispc
    mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -DNOCHECK -output ltv ltv.c
    mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -DNOCHECK -output ltv1Slice ltv1Slice.c
    mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -DNOCHECK -output cplxReassignment cplxReassignment.c
    mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -DNOCHECK -output fcwtComplexity fcwtComplexity.c
else
    disp('Platform not supported')
end