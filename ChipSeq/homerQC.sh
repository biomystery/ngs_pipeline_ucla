#!/bin/bash
while read data; do
    prefix=$(echo -e "$data"|cut -f1 -d".")
    makeTagDirectory "./homerQC/"$prefix -totalReads all -checkGC  $data &
done



