#!/bin/bash

PREFIX=marble.chkp

rm -f ${PREFIX}-??????.ppm

for chkpnt in `ls -1 ${PREFIX}-??????`
do
  printf 'Processing %s ... ' ${chkpnt}
  ./dump.x ${chkpnt} > ${chkpnt}.ppm
  printf 'Done\n'
done
