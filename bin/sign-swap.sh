#!/bin/bash
# ==========================================================================
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
# Authors: J-G. Dumas, B. Grenet, C. Pernet, A. Sedoglavic
# ==========================================================================
# ==========================================================================
# swaping the sign of a variable in an SLP
# ==========================================================================

var=$1
fil=$2
echo $0 ${var} ${fil}
sed "s/${var}/S${var}S/g;s/-S${var}S/M${var}/g;s/+S${var}S/P${var}/g;s/S${var}S/P${var}/g;s/M${var}/+${var}/g;s/P${var}/-${var}/g;s/:=+${var}/:=${var}/g;s/-${var}:=.*;/&)/;s/;)/);/;s/-${var}:=/${var}:=-(/;s/:=(.*)/PAPAPAPA&PAPAPAPA/;s/PAPAPAPA:=(/:=/;s/)PAPAPAPA//" ${fil} | compacter -s | sed 's/:=(.*)/PAPAPAPA&PAPAPAPA/;s/PAPAPAPA:=(/:=/;s/)PAPAPAPA//' | compacter -s | sponge ${fil}
