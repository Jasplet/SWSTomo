#! /bin/bash
############
# trig_bin.sh
############
# A program to take a given set of splitting results and geographically bin them
# using geogeom -T3. This is intended to be a lead up to using them splitting
# results with Matisse for SWS tomography
################################################################################

latlon=$1
echo "Binning Data at level T3"
echo $1
# cat $1 $2 |
