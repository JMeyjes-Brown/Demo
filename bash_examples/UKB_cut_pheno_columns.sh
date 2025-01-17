#! /bin/bash

COL_SELECT=$1
OUTPUT=$2

CUT_COLS=$(cut -f1$"{COL_SELECT}"\
             | grep -E [0-9]\
             | tr '\n' ',')

cut -f$"{CUT_COLS}" /mnt/work/datasets/compgen/ukbiobank/master/ukb43384.tab > $"{OUTPUT}"
