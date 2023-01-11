#!/bin/bash  -C

echo "tb2-ipr"

#toulbar2 bin path:
tb2="../build/bin/Linux"
# exec file : 
toulbar2="$tb2/toulbar2"
verbose=0


if [ -f "$toulbar2" ]; then

if [ "$TIMESEC" == "" ] ; then
	TIMESEC=9999999
fi

timer=$TIMESEC


if [ "$1" == "" ] ; then

if [ "$MODEL" == "" ] ; then
echo "$0 usage :  $0 foo.uai foo.uai.evid "
exit
fi
else 
MODEL=$1
fi 

if [ "$2" == "" ] ; then
	if [ "$EVID" == "" ] ; then
	EVID=$1.evid 
	echo "evid file set by default to : $EVID"
	fi

fi

if [ "$3" == "" ] ; then
	if [ "$RESULT" == "" ] ; then
	RESULT=$1.MPE
	echo "result set by default to : $RESULT"
	fi

fi


#echo "timesec $TIMESEC timer = $timer "

echo "$tb2/toulbar2  -timer=$timer  $OPT $MODEL $EVID"
#toulbar2 resolution 

else 
echo "$toulbar2 does not exist. please check tb2 path $tb2 "
exit
fi


# -e verbose=0 => toulbar2 output on screen
# -e verbose=1 => print debug info

dir=$tb2

#seed="-1"


if test -z "$verbose" 
then
      verbose="-1"
      echo "DEFAULT verbose LEVEL set to $verbose"
      
else
# verbose in  ENV
      echo "verbose LEVEL set to $verbose"
fi 

if test -z "$verbose" 
then
	verbose="-1"
fi

OPTIONE=" -i -pils -vns -O=-3 -E "
OPTION="  -A=100 -P=1000 -T=1000 -vacint -vacthr -rasps -raspsini -i -pils -p=-8 -O=-3 -t=1 -pwc=1 -n=10 -minqual "

SOLUTION=$RESULT.sol

if [ "$verbose" -ge "0" ] ; then
echo "timesec $TIMESEC";
echo "SOLUION FILE : $SOLUTION"
echo "MODEL=$MODEL"
echo "EVID=$EVID"
echo "RESULT=$RESULT"
echo "SOLUTION=$SOLUTION"

fi

if [ "$tsa" == "" ] ; then
	tsa=2
fi

if [ "$tsb" == "" ] ; then
	tsb=12
fi

if [ "$tsc" == "" ] ; then
	tsc=36
fi


timer=$TIMESEC

if [ "$TIMESEC" -eq 20 ] ; then
export	timer=$((TIMESEC-tsa)) 
fi

if [ "$TIMESEC" -eq 1200 ] ; then
export	timer=$((TIMESEC-tsb)) 
fi

if [ "$TIMESEC" -eq 3600 ] ; then
export	timer=$((TIMESEC-tsc)) 
fi


traj=$timer

# initial precision 2
# full precision : 7
pmax=2
for p in {4..7..3}
do
pmax=$((pmax+p))
done

#half time for HBFS-VPWC
pmax=$((pmax+pmax))

if [ "$verbose" -gt "0" ] ; then
echo "Number of intervals = $pmax"
fi

ut=0
#############
# initial solution with low precision 2
############
p=2
	

time=$(printf "%.0f\n" $(echo "scale=2; $((($p)*$TIMESEC/$pmax))" |bc))

if [ "$verbose" -gt "0" ] ; then
echo "---------------------"
echo "timer init $time"
echo "---------------------"
echo "used time $ut"
echo "---------------------"

echo "$dir/toulbar2 $MODEL $EVID -timer=$time  $OPTIONE -precision=2  -w=$RESULT -v=$verbose -seed=17101967 "
fi

if [ -f "$RESULT" ]; then
    rm -f "$RESULT"
if [ "$verbose" -ge "0" ] ; then
    echo "$RESULT deleted."
fi
fi

$dir/toulbar2 "$MODEL" "$EVID" -timer="$time" $OPTIONE -precision=2  -w="$RESULT" -v="$verbose" -seed=17101967

if [ -f "$SOLUTION" ]; then
    rm -f "$SOLUTION"
if [ "$verbose" -ge "0" ] ; then
    echo "$SOLUTION deleted."
fi
fi

awk 'FNR==2{$1="";print $0;exit}' "$RESULT" > "$SOLUTION"


if [ "$verbose" -ge "0" ] ; then
echo "SOLUTION $SOLUTION "
cat "$RESULT"
cat "$SOLUTION"
fi

tmp=$traj
traj=$(printf "%.0f\n" $(echo "scale=2; $(($tmp-$time))"|bc))
ut=$((ut+time))

for p in {4..7..3}
do

time=$(printf "%.0f\n" $(echo "scale=2; $((($p)*$TIMESEC/$pmax))" |bc))
tmp=$traj
traj=$(printf "%.0f\n" $(echo "scale=2; $(($tmp-$time))"|bc))
ut=$((ut+time))

if [ "$verbose" -ge "0" ] ; then
echo "------------"
echo " pres $p => duration $time"
echo "------------"

echo "$dir/toulbar2 $MODEL $EVID -timer=$time $OPTIONE -precision=$p $SOLUTION  -w=$RESULT -seed=$p"
fi

$dir/toulbar2 "$MODEL" "$EVID" -timer="$time" $OPTIONE -precision="$p" "$SOLUTION"  -w="$RESULT" -v="$verbose" -seed="$p"

if [ -f "$SOLUTION" ]; then
    rm -f "$SOLUTION"
#    echo "$SOLUTION deleted."
fi

awk 'FNR==2{$1="";print $0;exit}' "$RESULT" > "$SOLUTION"

if [ "$verbose" -ge "0" ] ; then
echo "SOLUTION:"
cat "$SOLUTION"
fi

done

res=$((timer-ut))


if [ "$verbose" -ge "0" ] ; then
echo "=========================="
echo " timer = $timer "
echo " used time $ut "
echo "proof search duration $res"
echo "=========================="

echo "$dir/toulbar2 -timer=$res $OPTION $MODEL $EVID $SOLUTION -w=$RESULT -v=$verbose "
fi

if [[ $res -gt 0 ]]; then
$dir/toulbar2 -timer=$res $OPTION "$MODEL" "$EVID" "$SOLUTION" -w="$RESULT" -v="$verbose"
fi

