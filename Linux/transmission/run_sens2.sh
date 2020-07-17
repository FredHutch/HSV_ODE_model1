#!/bin/sh
# This file launches a simulation multiple times until some number
# of transmissions occur.  A transmission is indicated by the following
# output message: "
#
# Parameters: <script> <dir> <input file> <# transmissions>
if [ "$1" != "" ]
then
if [ ! -d $1 ]
then
mkdir $1
fi
cd $1
else
echo "Please specify an input directory!"
fi

if [ "$2" != "" ]
then
cp ../$2 hsv_sim.in
cp ../hsv_sim.crit .
fi

if [ "$3" != "" ]
then
count=$3
else
count=500
fi

updates=`expr $count / 10`
next_update=`expr $count - $updates`

rm all_trans.out

while [ $count -gt 0 ]
do
echo "$count acts remaining"

../../hsv_sim  -c ../../hsv_sim.crit -r -b -w 2048 -v 2>&1 > trans.out
if [ ! -e "tcells.csv" ]
then
cp hsv_sim.dat12 tcells.csv
fi
trans=`grep -c "Transmission at" trans.out`
if [ $trans != "0" ]
then
acts=`grep "Transmission at" trans.out | cut -f2 -d'(' | cut -f1 -d' '`
count=`expr $count - $acts`
echo "Transmission after $acts acts" >> all_trans.out
else
echo "No transmissions after 500 acts" >> all_trans.out
count=`expr $count - 500`
fi
if [ $count -le $next_update ]
then
./send_update.sh
next_update=`expr $count - $updates`
fi
tail -125 trans.out >> all_trans.out
done
rm hsv_sim.dat12
