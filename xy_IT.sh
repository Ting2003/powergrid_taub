# generate output for xy[200-1000]
#!/bin/bash

BENCH_DIR=data/xy/wb
OUT_DIR=data/output/xy
WORKDIR=`pwd` 
PG=$WORKDIR/pg
MEMUSG=$WORKDIR/memusg
MODE=IT
PLATFORM=`uname`
TIME=/usr/bin/time
export LC_ALL=C

BENCHMARK=( 1000 800 600 500 400 300 250 200 )
#BENCHMARK=( 300 )
EPSILON=( 1e-5 ) #1e-4 1e-6 )
#OMEGA=( 1.4 1.5 1.6 1.7 1.8 1.9 1.3 1.2 1.1 1.0 ) #1.5 1.7 )
OMEGA=( 1.0 )
RATIO=( 0.2 )
BLOCKSIZE=( 7000 )
#TYPE=( y x )
TYPE=( y x )

# Start
for i in "${BENCHMARK[@]}"; do
for type in "${TYPE[@]}";do
NAME=$type$i.dat
INPUT=$WORKDIR/$BENCH_DIR/$NAME.spice.c4
BASENAME=`basename $INPUT`
for ep in "${EPSILON[@]}";do
for omega in "${OMEGA[@]}";do
for ratio in "${RATIO[@]}";do
for blocksize in "${BLOCKSIZE[@]}";do
	#RESULTDIR=$OUT_DIR/"$ep"_"$omega"_"$ratio"_"$blocksize"
	RESULTDIR=$OUT_DIR/
	mkdir -p $RESULTDIR
	OUTPUT=$WORKDIR/$OUT_DIR/$BASENAME.IT.out
	SYSTIME=$WORKDIR/$OUT_DIR/$BASENAME.IT.systime
	TIMELOG=$WORKDIR/$OUT_DIR/$BASENAME.IT.time

	# set options for time
	if [ "$PLATFORM" == 'Linux' ];then
		TIMEOPT="-v -p -o $SYSTIME"
	elif [ "$PLATFORM" == 'Darwin' ];then
		TIMEOPT="-l -p"
	fi

	#echo $NAME, $ep, $omega, $ratio, $blocksize
	echo $NAME, $omega
	
	CMD="$PG -e $ep -o $omega -r $ratio -b $blocksize -l $TIMELOG -i $INPUT -f $OUTPUT.tmp"
	#$TIME $TIMEOPT memusg $CMD $INPUT 2>$TIMELOG > $OUTPUT.tmp
	$TIME $TIMEOPT $MEMUSG $CMD 
	#memusg ls -alR / >/dev/null
	sort $OUTPUT.tmp > $OUTPUT
	rm $OUTPUT.tmp

	if [ -f $SYSTIME ];then
		cat $SYSTIME >> $TIMELOG
		rm $SYSTIME
	fi

	#mv -f "$TIMELOG" $RESULTDIR
	#mv -f "$OUTPUT" $RESULTDIR
done # blocksize
done # ratio
done # omega
done # ep
done # type
done # i
