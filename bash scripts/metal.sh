#!/usr/bin/env python

# Run make file for bins?
# Run Starlgiht
# turn starlight files into one fitted spectrum

#run twogauss
./twogauss < tginput

#run python
python metal.py
