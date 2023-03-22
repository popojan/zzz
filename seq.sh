#!/bin/bash

for i in `seq $1 $2`; do
  echo -n $i"	"
  ./zzz $i 0 0.1 10000 | tail -n1
done