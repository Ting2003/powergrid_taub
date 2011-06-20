#!/bin/bash

if [ "$#" -lt 1 ];then
	echo "Usage: ./batch.sh benchmark [mode]"
	exit 1
fi

export LC_ALL=C
ARGS="$@"

MODE=-I
if [ $MODE == '-L' ];then
	MODESTR=LU
else
	MODESTR=IT
fi

WORKDIR=`pwd` #/home/ac/tingyu1/powergrid
# Use Iterative method or LU
INPUT=$WORKDIR/$1
DIR=`dirname $INPUT`
BASENAME=`basename $INPUT`
OUTNAME=$DIR/$BASENAME.$MODESTR.out
SYSTIME=$DIR/$BASENAME.$MODESTR.systime
TIMELOG=$DIR/$BASENAME.$MODESTR.time

# set options for time
PLATFORM=`uname`
TIME=/usr/bin/time
if [ "$PLATFORM" == 'Linux' ];then
	TIMEOPT="-v -p -o $SYSTIME"
elif [ "$PLATFORM" == 'Darwin' ];then
	TIMEOPT="-l -p $TIMEOPT"
fi

MEMTOOL=$WORKDIR/memusg
PG=$WORKDIR/pg
# execute the program
$MEMTOOL $TIME $TIMEOPT $PG -l $TIMELOG $MODE -i $INPUT -f $OUTNAME

if [ -f $SYSTIME ];then
	cat $SYSTIME >> $TIMELOG
	rm $SYSTIME
fi

sort $OUTNAME > $OUTNAME.sort
mv $OUTNAME.sort $OUTNAME
