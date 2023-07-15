% minFunc
fprintf('Compiling minFunc files...\n');
if ismac || isunix
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -outdir minFunc/compiled minFunc/mex/mcholC.c
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -outdir minFunc/compiled minFunc/mex/lbfgsC.c
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS -ffp-model=fast" -v COPTIMFLAGS="-Ofast -fwrapv -DNDEBUG" -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c
elseif ispc
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -outdir minFunc/compiled minFunc/mex/mcholC.c
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -outdir minFunc/compiled minFunc/mex/lbfgsC.c
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
    mex -compatibleArrayDims -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c
else
    disp('Platform not supported')
end