#!/bin/bash  -C

echo "tb2-vns"
OPT=" -vns -O=-3 -E "


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

$toulbar2 -timer=$timer  $OPT $MODEL $EVID -w=$RESULT -v=$verbose
else 
echo "$toulbar2 does not exist. please check tb2 path $tb2 "
fi

