# sketchsort-minmax

Overview

## Description
SketchSort-minmax is a software for all pairs similarity search, a problem for 
finding all similar data pairs in a given dataset. 
We have implemented two different versions of SketchSort using different similarity measures; 
one uses Cosine similairy and the other uses Euclidean similarity, and 
both versions have been used for solving large-scale all pairs similarity search problems in various fields. 
This new version of SketchSort named SketchSort-minmax uses minmax as a similarity measure and is applicable to 
huge datasets as in the previous versions of SketchSort. 
SketchSort-minmax uses generalized consistent weighted sampling, a recently proposed random projection method for minmax similarity in KDD'17, 
for quickly solving large-scale all pairs similarity search. 

## Quick Start
cd src
make
./sketchsort-minmax -auto -znormalization -missingratio 0.0001 -minmax 0.1 ../dat/dat1000.txt outputfile

## Requirement
c++ compilar and boost library

## Usage
Usage: sketchsortj [OPTION]... INFILE OUTFILE

       where [OPTION]...  is a list of zero or more optional arguments

             INFILE       is the name of an input file

             OUTFILE      is the name of an output file

Additional arguments (input and output files may be specified):

       -hamdist [maximum hamming distance]

       (default: 1)

       -numblocks [the number of blocks]

       (default: 3)

       -minmax  [min-max distance threshold]

       (default: 0.1)

       -numchunks [the number of chunks]

       (default: 3)

       -auto 

       -missingratio [upper bound of false negative rate. large value makes sketchsort faster.]

       (default: 0.0001)

       -znormalization [dataset is z-nomalized]
       
       -minmaxnormalization [data is minmax-nomalized]

## Author
[Yasuo Tabei](https://sites.google.com/site/yasuotabei/)


