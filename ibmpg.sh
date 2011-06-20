# generate output for ibmpg[1-6]
#!/bin/bash

BENCH_DIR=data/ibmpg
OUT_DIR=data/output/ibmpg
WORKDIR=`pwd` 
PG=$WORKDIR/pg
MEMUSG=$WORKDIR/memusg
MODE=IT
PLATFORM=`uname`
TIME=/usr/bin/time
export LC_ALL=C

IBMPG=( 1 2 3 4 5 6 new1 new2 )
#IBMPG=( 6 )
#OMEGA=( 1.0 1.1 1.2 1.3 1.4 1.5 1.6 )
OMEGA=( 1.0 ) 

# Start
for i in "${IBMPG[@]}"; do
for omega in "${OMEGA[@]}";do
	NAME=ibmpg$i
	INPUT=$WORKDIR/$BENCH_DIR/$NAME.spice
	BASENAME=`basename $INPUT`
	#RESULTDIR=$WORKDIR/$OUT_DIR/omega_test
	RESULTDIR=$WORKDIR/$OUT_DIR
	SYSTIME=$RESULTDIR/$BASENAME.$MODE.systime
	TIMELOG=$RESULTDIR/$BASENAME.$MODE.time
	OUTPUT=$RESULTDIR/$BASENAME.$MODE.out
	#SYSTIME=$RESULTDIR/$BASENAME.$omega.systime
	#TIMELOG=$RESULTDIR/$BASENAME.$omega.time
	#OUTPUT=$RESULTDIR/$BASENAME.$omega.out
	mkdir -p $RESULTDIR
	# set options for time
	if [ "$PLATFORM" == 'Linux' ];then
		TIMEOPT="-v -p -o $SYSTIME"
	elif [ "$PLATFORM" == 'Darwin' ];then
		TIMEOPT="-l -p"
	fi

	echo $NAME
	
	CMD="$PG -l $TIMELOG -i $INPUT -f $OUTPUT.tmp"
	#$TIME $TIMEOPT memusg $CMD $INPUT 2>$TIMELOG > $OUTPUT.tmp
	$TIME $TIMEOPT $MEMUSG $CMD
	sort $OUTPUT.tmp > $OUTPUT
	rm $OUTPUT.tmp

	if [ -f $SYSTIME ];then
		cat $SYSTIME >> $TIMELOG
		rm $SYSTIME
	fi
done #OMEGA
done # i
