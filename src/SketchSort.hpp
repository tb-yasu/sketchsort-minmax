/*
 * SketchSort.hpp
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

#ifndef SKETCHSORT_HPP
#define SKETCHSORT_HPP

// for boost
#include <boost/pool/pool.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/random.hpp>

#include <random>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <cstring>
#include <string>
#include <iterator>
#include <fstream>
#include <strstream>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <stdint.h>
#include <time.h>

#ifndef DBL_MAX
  #define  DBL_MAX 1.797693e+300
#endif

class Random{
public:
  Random() {
    srand( static_cast<uint32_t>( time(NULL) ) );
  }
  uint32_t operator()(uint32_t max) {
    double tmp = static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
    return static_cast<uint32_t>( tmp * max );
  }
};


struct pstat {
  int start;
  int end;
};

struct params {
  uint32_t numblocks;
  uint32_t numchunks;
  uint32_t projectDim;
  uint32_t chunk_dist;
  uint32_t chunks;
  uint32_t num_seq;
  uint32_t seq_len;
  uint32_t chunk_len;
  uint32_t start_chunk;
  uint32_t end_chunk;
  uint32_t cchunk;
  uint32_t *counter;
  pstat        *pos;
  pstat        *pchunks;
  float         minmaxDist;
  std::vector<uint32_t> ids;
  std::vector<int> blocks;
  std::ostream *os;
};

class SketchSort {
private:
  void readFeature(const char *fname);
  void projectVectors(uint32_t projectDim, std::vector<uint8_t*> &sig, params &param);
  void report(std::vector<uint8_t*> &sig, int l, int r, params &param);
  void radixsort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param);
  void insertionSort(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param);
  bool calc_chunk_hamdist(uint8_t *seq1, uint8_t *seq2, const params &param);
  bool calc_hamdist(uint8_t *seq1, uint8_t *seq2, const params &param);
  bool check_canonical(uint8_t *seq1, uint8_t *seq2, const params &param);
  bool check_chunk_canonical(uint8_t *seq1, uint8_t *seq2, const params &param);
  float calcMinMaxDist(uint32_t id1, uint32_t id2);
  float check_upper_bound(uint32_t id1, uint32_t id2);
  double calcMissingEdgeRatio(params &param);
  void classify(std::vector<uint8_t*> &sig, int spos, int epos, int l, int r, int bpos, params &param);
  void turnParameters(float _missingratio, params &param);
  void multi_classification(std::vector<uint8_t*> &sig, int maxind, int l, int r, params &param);

  float log10(float val);
  void generateMatrix(uint32_t projectDim, std::vector<std::vector<float> > &matR, std::vector<std::vector<float> > &matC, std::vector<std::vector<float> > &matB);
  void convertVector(std::vector<float> &fv, std::vector<float> &fv2);
  uint64_t compMin(std::vector<float> &resA);
  void zNormalization();
  void MinMaxNormalization();
public:
  SketchSort() {};
  void run(const char *fname, const char *oname, 
	   uint32_t _numblocks, 
	   uint32_t _dist, 
	   float        _minmaxDist,
	   uint32_t _numchunks, 
	   bool         _autoFlag,
	   float        _missingratio, 
	   bool         _zNormalization, 
	   bool         _minmaxNormalization); 
private:
  boost::pool<> *p;
  std::vector<std::vector<float> > fvs;
  uint64_t           dim_;
  uint32_t           num_char;
};

#endif // SKETCHSORT_HPP
