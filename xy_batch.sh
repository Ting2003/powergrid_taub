# generate output for xy[200-1000]
#!/bin/bash

BENCH_DIR=data/xy
OUT_DIR=data/output/xy/LU

for i in 1000 800 600 500 400 300 250 200; 
do
	XNAME=x$i.dat
	YNAME=y$i.dat
	X=$BENCH_DIR/$XNAME.spice
	Y=$BENCH_DIR/$YNAME.spice
	X_OUT=$OUT_DIR/$XNAME.out
	Y_OUT=$OUT_DIR/$YNAME.out

	echo "processing $X ..."
	./batch.sh $X
	#if [ "$?" -eq 1 ];then 
	#	echo "Fail to continue"
	#	exit 1
	#fi

	echo "processing $Y ..."
	./batch.sh $Y
done
