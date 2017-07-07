#!/bin/bash

PREFIX=marble.chkp

rm -f ${PREFIX}-??????.ppm ${PREFIX}-??????.png

for chkpnt in `ls -1 ${PREFIX}-??????`
do
  printf 'Processing %s ... ' ${chkpnt}
  ./dump.x ${chkpnt} > ${chkpnt}.ppm
  convert ${chkpnt}.ppm -resize 150x150 ${chkpnt}.png
  printf 'Done\n'
done
