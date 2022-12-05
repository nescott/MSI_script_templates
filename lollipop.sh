#! /bin/bash

gene=("$@")
~/bin/lollipops-v1.5.2-linux64/lollipops -labels -U ${gene[0]} -l=${gene[1]} ${gene[@]:2:($#-1)}

echo "$#"
echo "$@"
count=($# - 1)
echo $count
