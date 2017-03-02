#!/usr/bin/env python
#!/bin/bash

# Run make file for bins?
massbin = $1
redbin = $2

sh ./spectra/wig02make.txt $massbin $redbin

# Run Starlgiht
./home/uqzbyrne/Documents/STARLIGHTv04/StarlightChains_v04.exe < wig024363.in > wig02o3w.out
./home/uqzbyrne/Documents/STARLIGHTv04/StarlightChains_v04.exe < wig02o2.in > wig02o2.out
./home/uqzbyrne/Documents/STARLIGHTv04/StarlightChains_v04.exe < wig02hb.in > wig02hb.out


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
