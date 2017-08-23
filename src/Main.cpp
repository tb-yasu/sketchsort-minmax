/*
 * Main.cpp
 * Copyright (c) 2017 Yasuo Tabei All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE and * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SketchSort.hpp"

#include <iostream>
#include <cstdlib>

using namespace std;

/* Globals */
void usage();
void version();
void parse_parameters (int argc, char **argv);

char *fname, *oname;
int   numblocks           = 3;
int   hamDist             = 1;
int   numchunks           = 3;
float minmaxDist          = 0.1;
bool  autoFlag            = false;
float missingratio        = 0.0001;
bool  zNormalization      = false;
bool  minmaxNormalization = false;

int main(int argc, char **argv) 
{
  version();

  parse_parameters(argc, argv);

  SketchSort sketchsort;
  sketchsort.run(fname, oname, numblocks, hamDist, minmaxDist, numchunks, autoFlag, missingratio, zNormalization, minmaxNormalization);

  return 0;
}

void version(){
  cerr << "SketchSort Min-Max Distance version 0.0.1" << endl;
  cerr << "Written by Yasuo Tabei" << endl << endl;
}

void usage(){
  cerr << endl
       << "Usage: sketchsortj [OPTION]... INFILE OUTFILE" << endl << endl
       << "       where [OPTION]...  is a list of zero or more optional arguments" << endl
       << "             INFILE       is the name of an input file" << endl
       << "             OUTFILE      is the name of an output file" << endl << endl
       << "Additional arguments (input and output files may be specified):" << endl
       << "       -hamdist [maximum hamming distance]" << endl
       << "       (default: " << hamDist << ")" << endl
       << "       -numblocks [the number of blocks]" << endl
       << "       (default: " << numblocks << ")" << endl
       << "       -minmax  [min-max distance threshold]" << endl
       << "       (default: " << minmaxDist << ")" << endl
       << "       -numchunks [the number of chunks]" << endl
       << "       (default: " << numchunks << ")" << endl
       << "       -auto " << endl
       << "       -missingratio" << endl
       << "       -znormalization" << endl
       << "       -minmaxnormalization"  << endl
       << "       (default: " << missingratio << ")" << endl
       << endl;
  exit(0);
}

void parse_parameters (int argc, char **argv){
  //  if (argc == 1) usage();
  int argno;
  for (argno = 1; argno < argc; argno++){
    if (argv[argno][0] == '-'){
      if      (!strcmp (argv[argno], "-version")){
	version();
      }
      else if (!strcmp (argv[argno], "-numblocks")) {
	if (argno == argc - 1) cerr << "Must specify minimum support after -numblocks" << endl;
	numblocks = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-auto")) {
	autoFlag = true;
      }
      else if (!strcmp (argv[argno], "-hamdist")) {
	if (argno == argc - 1) cerr << "Must specify miximum itemset size after -hamdist" << endl;
	hamDist = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-minmax")) {
	if (argno == argc - 1) cerr << "Must specify miximum itemset size after -minmax" << endl;
	minmaxDist = atof(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-numchunks")) {
	if (argno == argc - 1) cerr << "Must specify minimum support after -numchunks" << endl;
	numchunks = atoi(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-missingratio")) {
	if (argno == argc - 1) cerr << "Must specify missing edge ratio after -missingratio" << endl;
	missingratio = atof(argv[++argno]);
      }
      else if (!strcmp (argv[argno], "-znormalization")) {
	zNormalization = true;
      }
      else if (!strcmp (argv[argno], "-minmaxnormalization")) {
	minmaxNormalization = true;
      }
      else {
	usage();
      }
    } else {
      break;
    }
  }
  if (argno + 1 >= argc)
    usage();

  fname = argv[argno];
  oname = argv[argno + 1];
}
