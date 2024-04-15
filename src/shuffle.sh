#!/bin/bash 

cd $1

fs=(*/str.txt)
sfs=($(shuf -e "${fs[@]}"))

size=${#fs[*]}

for (( i=0; $i < $size; i++ ))
do
	mv ${fs[$i]} ${sfs[$i]}.new
done

for f in ./*/str.txt.new
do 
	mv -- "$f" "${f%.new}"
done 

