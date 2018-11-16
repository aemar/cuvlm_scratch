#!/bin/bash
#make clean
#make vlatu
for i in {1..5};
do
cp -v -r input/job$i/{vlatinit.dat,vlatin0001.dat,trajectory0001.dat} input
./vlatu
mkdir input/job$i/output
cp -v -r output/{cp,wake,wing,cl.dat,cd.dat,rundata.dat,wake.ply,wing.ply,wing_trnsfrm.ply} input/job$i/output
#(cd input/job$i; ../../transform/fftcp)
#if [ $i -eq 4 ]; then
#make clean
#make vlatu2
#fi
done
