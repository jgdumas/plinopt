#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# Script: finds 3 files with starting letters L, R and P
# ==========================================================================

basename $1 | awk -v dir=`dirname $1` '{ c=substr($1, 2); print dir"/L"c,dir"/R"c,dir"/P"c }'
