#!/bin/sh

# Credit: 
#       1. Zheng,W., Wuyun,Q., et al. (2019) Detecting distant-homology protein structures by aligning deep neural-network based contact maps. PLOS Computational Biology

alnfile=$1
a3mfile=$2
i=0;
while read line
do
if [ "$i" -eq "0" ]
then
   echo ">protein" >$a3mfile
   echo "$line"    >>$a3mfile
   i=`expr $i + 1`
else
   echo ">seq_$i"  >>$a3mfile
   echo "$line"   >>$a3mfile
   i=`expr $i + 1`
fi
done<$alnfile
