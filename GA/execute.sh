#!/bin/bash
for filename in instances/teste/*.tsp; do
    for ((i=0; i<10; i++)); do
        (printf %s "$(cut -d'/' -f3 <<< $filename)"; printf ", "; ./tsp "$filename" ) | cat
    done
done