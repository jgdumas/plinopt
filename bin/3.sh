#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# Script: finds 3 files ending with letters L, R and P then .sms
# ==========================================================================

basename $1 .sms | awk -v dir=`dirname $1` '{ c=substr($1,1,length($1)-1); print dir"/"c"L.sms",dir"/"c"R.sms",dir"/"c"P.sms" }'
