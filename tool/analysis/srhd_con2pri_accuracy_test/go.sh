#!/bin/bash

NAME=$1


if [ -z "$NAME" ];then
   echo "error: Please provide source code"
   exit
fi


OUT=${NAME%.cpp}

PRECISION=FLOAT8


g++ $1  CPU_Shared_FluUtility.cpp \
CPU_EoS_TaubMathews.cpp \
EoS_Init.cpp \
Aux_Error.cpp \
Aux_Message.cpp \
-DSERIAL \
-DLR_SCHEME=PLM \
-DMODEL=HYDRO \
-DSRHD \
-DRANDOM_NUMBER=RNG_GNU_EXT \
-DFLU_SCHEME=MHM \
-DRSOLVER=HLLC \
-DNLEVEL=10 \
-DMAX_PATCH=200000 \
-DREDUCED_ENERGY \
-D$PRECISION \
-DEOS=EOS_TAUBMATHEWS -lm -Wall -o $OUT && ./$OUT
