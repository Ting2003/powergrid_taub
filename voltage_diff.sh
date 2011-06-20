#!/bin/bash

if [ "$#" -lt 2 ];then
	echo "Usage: ./voltage_diff.sh file1 file2"
	exit 1
fi

FILE1=$1
FILE2=$2
#FILE1=data/output/ibmpg/ibmpg3.out
#FILE2=data/ibmpg/solution/ibmpg3.solution.sort 
n1=`wc -l $FILE1 | cut -d ' ' -f 1`
n2=`wc -l $FILE2 | cut -d ' ' -f 1`
if [[ $n1 != $n2 ]];then
	echo "Two files have different line numbers!"
	exit 1
fi
join $FILE1 $FILE2 | awk -f voltage_diff.awk
