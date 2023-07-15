mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output initCplxMovAvg initCplxMovAvg.c
mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output procCplxMovAvg procCplxMovAvg.c
mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output initNPole1Zero initNPole1Zero.c
mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output procNPole1Zero procNPole1Zero.c
mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output freeFDL freeFDL.c
mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output initFDL initFDL.c
mex -R2018a -v COMPFLAGS="$COMPFLAGS /GL /fp:fast" -output procFDL procFDL.c