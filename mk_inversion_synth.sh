#! /bin/bash

# Parse arguement
# expected usage mk_inversion_synth.sh -spol 30 -n 0.01 -file filename.sdb
POS=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
  -s|-spol)
  SPOL=$2
  shift # past arguement (shift makes arg $2 become arg $1 etc (shifts args along 1))
  shift # past value
  ;;
  -n|-noise)
  NOISE_LVL=$2
  shift
  shift
  ;;
  -f|-file)
  SDB=$2
  shift # past -file flag
  shift # past file
  ;;
  *)
  POS+=("$1") # Any extra args
  shift
  ;;
esac
done

function call_sacsplitwave {
#    Function to basically call sacsplitwav
echo $1 $2 $3 $4 $5 $6 #var $5 and $6 are for upper layer if provided
if [ -z "$5" ]
then
  sacsplitwave -op $1 $2 -spol $3 -dfreq 0.125 -noise $4
else
  sacsplitwave -op $1 $2 -op $5 $6 -spol $3 -dfreq 0.125 -noise $4
fi
#Use 0.1 for "low noise", 0.25 for "high" noise and now also add "0.05" for "very low"
}

i=1 
while read file; do
    echo $i
    #SDB *should* have headers 
    #STAT DATE TIME PHASE EVLA EVLO EVDP STLA STLO DIST BAZ AZI INCL LENGTH FAST TLAG
    echo $file
    if [ $i -gt 1 ]
    then
	    call_sacsplitwave $(echo $file | awk '{print $18, $19}' ) $SPOL `echo "scale=2;$NOISE_LVL/100" | bc `
 	stat=`echo $file | awk '{print $1}'`
 	date=`echo $file | awk '{print $2}'`
 	time=`echo $file | awk '{print $3}'`
  ph=`echo $file | awk '{print $4}'`
 	path=`pwd`

  for cmp in E N Z
  do
   	mv SWAV.BH${cmp} ${path}/data/${stat}_synth_${ph}.BH${cmp}
    # Alter station name, Date
    sacsethdr kstnm ${stat} ${path}/data/${stat}_synth_${ph}.BH${cmp}
    sacsethdr nzyear `echo ${date} | cut -c 1-4` ${stat}_synth_${ph}.BH${cmp}
    sacsethdr kcmpnm BH${cmp} ${stat}_synth_${ph}.BH${cmp}
  # sacsethdr nztime 12${snr} SP${SPOL}_${k}001_120000.BH[E,N,Z]
  echo "${path}/${stat}_synth_${ph}" | cat >> Inversion_synthetics_noise_${NOISE_LVL}.events
	else
  echo 'Skipping Header Line'
  fi
  let i=i+1
done < $SDB
