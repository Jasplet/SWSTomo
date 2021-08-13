#!/usr/bin/env bash
######################
# bin_by_line.sh
######################
# A short script to use geogeom to
# bin ScS line by line, so that we can find out
# which line geogeom does not like
path=/Users/ja17375/DiscrePy/Sheba/Results/ScS
outfile=$1
# rm $path/bins.out
# cat $path/bins.out
# Make Header line
echo "BIN MID_LAT MID_LON V1_LAT V1_LON V2_LAT V2_LON V3_LAT V3_LON" > $outfile
i=0
while read line
do
  if [ $i -gt 0 ] # i.e skip header line
  then
  # echo $i >> $path/bins.out
  echo $line | awk '{print $4,$5}' | geogeom T3+ -- > tmp
  printf "%d %f %f %f %f %f %f %f %f\n" `cat tmp` >> $outfile
  fi
  let i+=1
done < $path/ScS_w_bouncepoints.sdb
rm tmp
