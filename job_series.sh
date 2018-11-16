#!/bin/bash
#make clean
#make vlatu_notree
#for i in {12}; do
  for j in {05,1,2,3,4,5}; do
    (cd input/cl_M_var_new/AR12/M0$j; ../../../../vlatu
    cp output/cl.dat ../../cl.AR12.M0$j.dat) &
  done
#done

#05,1,2,3,4,5,6,7,8,85,9,95
