# generate output for xy[200-1000]
#!/bin/bash

export LC_ALL=C
BENCH_DIR=data/xy
OUT_DIR=data/output/xy
WORKDIR=`pwd` 
PG=$WORKDIR/IPGSv1.1
MEMUSG=$WORKDIR/memusg
MODE=LU
PLATFORM=`uname`
TIME=/usr/bin/time

#BENCHMARK=( 1000 800 600 500 400 300 250 200 )
BENCHMARK=( 200 250 )
TYPE=( y x )

# Start
for i in "${BENCHMARK[@]}"; do
for type in "${TYPE[@]}";do
	NAME=$type$i.dat
	INPUT=$WORKDIR/$BENCH_DIR/$NAME.spice
	BASENAME=`basename $INPUT`
	OUTPUT=$WORKDIR/$OUT_DIR/$BASENAME.$MODE.out
	SYSTIME=$WORKDIR/$OUT_DIR/$BASENAME.$MODE.systime
	TIMELOG=$WORKDIR/$OUT_DIR/$BASENAME.$MODE.time
	RESULTDIR=$OUT_DIR
	mkdir -p $RESULTDIR
	# set options for time
	if [ "$PLATFORM" == 'Linux' ];then
		TIMEOPT="-v -p -o $SYSTIME"
	elif [ "$PLATFORM" == 'Darwin' ];then
		TIMEOPT="-l -p"
	fi

	echo $NAME
	
	CMD="$PG -L -l $TIMELOG -i $INPUT -f $OUTPUT.tmp"
		
	$TIME $TIMEOPT $MEMUSG $CMD
	sort $OUTPUT.tmp > $OUTPUT
	rm $OUTPUT.tmp

	if [ -f $SYSTIME ];then
		cat $SYSTIME >> $TIMELOG
		rm $SYSTIME
	fi

done # type
done # i
