#!/bin/bash
n=$1
t=0
all_t=""

mkdir -p out/

printf "%12s\t%12s\n" "seed" "user-time (s)"

TIMEFORMAT=%U
for i in $(seq 0 9); do
    t=$( { time bin/hw2.run --gtsp -s $i -n $n 2>&1 > out/gtsp-n${n}-s${i}.log; } 2>&1 )
    all_t="$all_t $t"
    printf "%12d\t%12.3f\n" $i $t
done

echo
echo $all_t | scripts/mean.r
