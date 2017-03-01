#!/usr/bin/env python
#!/bin/bash

# Run make file for bins?
massbin = $1
redbin = $2

sh ./wig02make.txt $massbin $redbin

# Run Starlgiht


# turn starlight files into one fitted spectrum
for i in {0..massbin}
do
for j in {0..redbin}
do
python create2gauss.py $i$j
done
done

ls wig02fitted*.txt > list.txt

#run twogauss
./twogauss < tginput

#run python
python metal.py
