# sketchsort-minmax

Overview

## Description
SketchSort-minmax is a software for all pairs similarity search, a problem for 
finding all similar data pairs from a given dataset. 
We have implemented two different versions of SketchSort using different similarity measures; 
one uses Cosine similairy and the other uses Euclidean similarity, and 
SketchSort has been used for solving large-scale all pairs similarity search problems in various fields. 
This new version of SketchSort named SketchSort-minmax uses minmax as a similarity measure and is applicable to 
huge datasets as in the previous versions of SketchSort. 
SketchSort-minmax uses generalized consistent weighted sampling, a recently proposed random projection method for minmax similarity in KDD'17, 
for quickly solving large-scale all pairs similarity search. 

## Quick Start
cd src
make


