/*
 * SketchSort.cpp
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

using namespace std;
using namespace boost;

float SketchSort::log10(float val) {
  if (val < 1.0e-20) {
    return -1.0e+20;
  }

  return log(val);
}

void SketchSort::readFeature(const char *fname) {
  ifstream ifs(fname);

  if (!ifs) {
    std::cerr << "can not open " << fname << std::endl;
    exit(0);
  }

  float val = 0;
  string line;
  while (std::getline(ifs, line)) {
    if (line.size() == 0)
      continue;
    fvs.resize(fvs.size() + 1);
    vector<float> &fv = fvs[fvs.size() - 1];
    istringstream is(line);
    if (fvs.size() == 1) {
      while (is >> val)
	fv.push_back(val);
      fv.swap(fv);
      dim_ = fv.size();
    }
    else {
      fv.resize(dim_);
      uint32_t iter = 0;
      while (is >> val) {
	if (iter >= dim_) {
	  cerr << "error: data dimension is different." << endl;
	  exit(1);
	}
	fv[iter++] = val;
      }
    }
  }
}

void SketchSort::generateMatrix(uint32_t projectDim, vector<vector<float> > &matR, vector<vector<float> > &matC, vector<vector<float> > &matB) {
  uint64_t datDim2 = dim_ << 1;
  
  {
    boost::mt19937 gen1(static_cast<unsigned long>(time(0)));
    boost::gamma_distribution<> dst1(2.0);
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > grand1(gen1, dst1);

    matR.resize(projectDim);
    for (size_t i = 0; i < projectDim; ++i) {
      vector<float> &vecR = matR[i];
      vecR.resize(datDim2);
      for (size_t j = 0; j < datDim2; ++j) {
	vecR[j] = grand1();
      }
    }
  }
  
  {
    boost::mt19937 gen2(static_cast<unsigned long>(time(0)));
    boost::gamma_distribution<> dst2(2.0);
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > grand2(gen2, dst2);

    matC.resize(projectDim);
    for (size_t i = 0; i < projectDim; ++i) {
      std::vector<float> &vecC = matC[i];
      vecC.resize(datDim2);
      for (size_t j = 0; j < datDim2; ++j) 
	vecC[j] = grand2();
    }
  }

  {
    boost::mt19937        gen3(static_cast<uint32_t>(time(0)));
    boost::uniform_real<> dst3(0, 1.0);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > urand(gen3, dst3);
    
    matB.resize(projectDim);
    for (size_t i = 0; i < projectDim; ++i) {
      vector<float> &vecB = matB[i];
      vecB.resize(datDim2);
      for (size_t j = 0; j < datDim2; ++j) {
	vecB[j] = urand();
      }
    }
  }
}

void SketchSort::convertVector(vector<float> &fv, vector<float> &fv2) {
  fv2.resize(dim_ << 1);

  uint64_t iter = 0;
  for (size_t i = 0; i < dim_; ++i) {
    if (fv[i] > 0) {
      fv2[iter++] = fv[i];
      fv2[iter++] = 0.0;
    }
    else {
      fv2[iter++] = 0.0;
      fv2[iter++] = -fv[i];
    }
  }
}

uint64_t SketchSort::compMin(vector<float> &resA) {
  float minVal = 3.402823466e+35F;
  uint64_t index = 0;
  for (size_t i = 0; i < resA.size(); ++i) {
    if (resA[i] < minVal) {
      minVal = resA[i];
      index = i;
    }
  }
  
  return index;
}

void SketchSort::projectVectors(uint32_t projectDim, vector<uint8_t*> &sig, params &param) {
  vector<vector<float> > matR;
  vector<vector<float> > matC;
  vector<vector<float> > matB;
  generateMatrix(projectDim, matR, matC, matB);

  p = new boost::pool<>(sizeof(uint8_t));
  sig.resize(fvs.size());
  param.ids.resize(fvs.size());
  for (size_t i = 0; i < sig.size(); i++) {
    sig[i]       = (uint8_t*)p->ordered_malloc(projectDim + 1);
    param.ids[i] = i;
  }

  uint64_t datDim2 = dim_ << 1;

  for (size_t i = 0; i < fvs.size(); ++i) {
    if (i % 1000 == 0)
      cerr << "project " << i << " < " << fvs.size() << endl;
    
    vector<float> &fv = fvs[i];

    for (size_t j = 0; j < projectDim; ++j) {
      vector<float> &vecR = matR[j];
      vector<float> &vecC = matC[j];
      vector<float> &vecB = matB[j];

      float minVal = 3.402823466e+35F;
      uint64_t minIndex = 0;
      for (size_t k = 0; k < dim_; ++k) {
	float val1 = 0.0;
	float val2 = 0.0;
	
	if (fv[k] > 0) 
	  val1 = fv[k];
	if (fv[k] < 0)
	  val2 = -fv[k];

	uint64_t ind = (k << 1);
	int64_t t = floor(log10(val1)/vecR[ind] + vecB[ind]);
	float a = log10(vecC[ind]) - vecR[ind] * ((float)t + 1.0 - vecB[ind]);
	if (a < minVal) {
	  minVal = a;
	  minIndex = ind;
	}

	ind = (k << 1) + 1;
	t = floor(log10(val2)/vecR[ind] + vecB[ind]);
	a = log10(vecC[ind]) - vecR[ind] * ((float)t + 1.0 - vecB[ind]);
	if (a < minVal) {
	  minVal = a;
	  minIndex = ind;
	}
      }
      sig[i][j+1] = minIndex%256;
    }
  }
  
  param.seq_len = projectDim;
  param.num_seq = fvs.size();

}

inline float SketchSort::calcMinMaxDist(uint32_t id1, uint32_t id2) {
  vector<float> &fv1 = fvs[id1];
  vector<float> &fv2 = fvs[id2];

  double num = 0.0, den = 0.0;
  for (size_t d = 0; d < dim_; ++d) {
    if (fv1[d] > 0 && fv2[d] > 0) {
      num += min(fv1[d], fv2[d]);
      den += max(fv1[d], fv2[d]);
    }
    else if (fv1[d] < 0 && fv2[d] < 0) {
      num += min(-fv1[d], -fv2[d]);
      den += max(-fv1[d], -fv2[d]);
    }
    else {
      den += fabs(fv1[d]);
      den += fabs(fv2[d]);
    }
  }

  return (1.f - num/den);
}

inline void SketchSort::radixsort(vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param) {
  uint32_t *c      = param.counter;
  vector<uint32_t> &ids      = param.ids;
  vector<uint8_t*> newsig(r - l + 1);
  vector<uint32_t> newids(r - l + 1);
  uint32_t tmp;
  int tpos = spos - 1;
  while (++tpos <= epos) {
    for (int i = 0; i < num_char; i++) *(c + i) = 0;
    for (int i = l; i <= r; i++) c[sig[i][tpos]]++;
    for (int i = 1; i < num_char; i++) *(c + i) += *(c + i - 1);
    for (int i = r; i >= l; --i) {
      tmp = --c[sig[i][tpos]] + l;
      newids[tmp - l] = ids[i];
      newsig[tmp - l] = sig[i];
    }
    if (++tpos <= epos) {
      for (int i = 0; i < num_char; i++) *(c + i) = 0;
      for (int i = l; i <= r; i++) c[newsig[i - l][tpos]]++;
      for (int i = 1; i < num_char; i++) *(c + i) += *(c + i - 1);
      for (int i = r; i >= l; --i) {
	tmp = --c[newsig[i - l][tpos]] + l;
	ids[tmp] = newids[i - l];
	sig[tmp] = newsig[i - l];
      }
    }
    else {
      for (int i = l; i <= r; i++) {
	ids[i] = newids[i - l];
	sig[i] = newsig[i - l];
      }
      return;
    }
  }
}

inline void SketchSort::insertionSort(vector<uint8_t*> &sig, int spos, int epos, int l, int r, params &param) {
  int i, j;
  uint8_t *pivot, pval;
  uint32_t pid;
  vector<uint32_t> &ids = param.ids;
  for (int tpos = spos; tpos <= epos; tpos++) {
    for (i = l + 1; i <= r; i++) {
      pivot = sig[i]; pval = sig[i][tpos]; pid = ids[i];
      for (j = i; j > l && sig[j-1][tpos] > pval; j--) {
	sig[j] = sig[j-1];
	ids[j] = ids[j-1];
      }
      sig[j] = pivot;
      ids[j] = pid;
    }
  }
}

inline void SketchSort::classify(vector<uint8_t*> &sig, int spos, int epos, int l, int r, int bpos, params &param) {
  int n_l = l, n_r = r;
  for (int iter = l + 1; iter <= r; iter++) {
    if (!equal(sig[n_l] + spos, sig[n_l] + epos + 1, sig[iter] + spos)) {
      n_r = iter - 1;
      if (n_r - n_l >= 1)
	multi_classification(sig, bpos + 1, n_l, n_r, param);
      n_l = iter;
    }
  }
  if (r - n_l >= 1)
    multi_classification(sig, bpos + 1, n_l, r, param);
}

inline bool SketchSort::calc_chunk_hamdist(uint8_t *seq1, uint8_t *seq2, const params &param) {
  uint32_t d = 0;
  for (size_t i = 1;  i <= param.chunk_len; i++) 
    if (*seq1++ != *seq2++ && ++d > param.chunk_dist) return false;
  return true;
}

inline bool SketchSort::check_chunk_canonical(uint8_t *seq1, uint8_t *seq2, const params &param) {
  uint32_t d = 0;
  int end        = param.pchunks[param.cchunk].start - 1;
  int j          = 1;
  int tend       = param.pchunks[j].end;
  int i          = 0;
  
  while (++i <= end) {
    if (seq1[i] != seq2[i] && ++d > param.chunk_dist) {
      while (++i <= tend)
	if (seq1[i] != seq2[i]) ++d;
      d = 0;
      tend = param.pchunks[++j].end;
      i    = param.pchunks[j].start - 1;
      continue;
    }
    if (tend == i)
      return false;
  }
  return true;
}

inline bool SketchSort::check_canonical(uint8_t *seq1, uint8_t *seq2, const params &param) {
  size_t sb = 1, eb = 1;
  size_t b;
  for (size_t i = 0, size = param.blocks.size(); i < size; i++) {
    eb = param.blocks[i];
    for (b = sb; b < eb; b++) {
      if (equal(seq1 + param.pos[b].start, seq1 + param.pos[b].end + 1, seq2 + param.pos[b].start))
	  return false;
    }
    sb = param.blocks[i] + 1;
  }
  return true;
}


inline void SketchSort::report(vector<uint8_t*> &sig, int l, int r, params &param) {
  //  std::cout << "report" << std::endl;
  float minmaxDist;
  for (int i = l; i < r; i++) {
    for (int j = i + 1; j <= r; j++) {
      if (calc_chunk_hamdist(sig[i] + param.start_chunk, sig[j] + param.start_chunk, param) && 
	  check_canonical(sig[i], sig[j], param) && 
	  check_chunk_canonical(sig[i], sig[j], param) &&
	  ((minmaxDist = calcMinMaxDist(param.ids[i], param.ids[j])) <= param.minmaxDist)) {
	(*param.os) << param.ids[i] << " " << param.ids[j] << " " << minmaxDist << endl;
      }
    }
  }
}

void SketchSort::multi_classification(vector<uint8_t*> &sig, int maxind, int l, int r, params &param) {
  if (param.blocks.size() == param.numblocks - param.chunk_dist) {
    report(sig, l, r, param);
    return;
  }

  for (int bpos = maxind; bpos <= (int)param.numblocks; bpos++) {
    if (param.blocks.size() + (param.numblocks - bpos + 1) < param.numblocks - param.chunk_dist) { // pruning
      //      std::cerr << "return " << std::endl;
      return;
    }
    param.blocks.push_back(bpos);
    radixsort(sig, param.pos[bpos].start, param.pos[bpos].end, l, r, param);
    classify(sig, param.pos[bpos].start, param.pos[bpos].end, l, r, bpos, param);
    param.blocks.pop_back();
  }
}

double combination(int n, int m) {
  double sum = 1.0;
  for (int i = 0; i < m; i++) {
    sum *= ((n-i)/(m-i));
  }
  return sum;
}

double SketchSort::calcMissingEdgeRatio(params &param) {
  double sum = 0.f;
  double prob = (255.0/256.0)*param.minmaxDist;
  uint32_t l = param.projectDim;
  uint32_t d = param.chunk_dist;
  for (uint32_t k = 0; k <= d; k++) {
    sum += (combination(l, k) * pow(prob, k) * pow(1.f - prob, l - k));
  }
  return pow(1.f - sum, param.numchunks);
}

void SketchSort::turnParameters(float _missingratio, params &param) {
  uint32_t hamDist   = 1;
  uint32_t numBlocks = hamDist + 3;
  uint32_t numchunks = 0;

  do {
    if (numchunks > 30) {
      hamDist   += 1;
      numBlocks  = hamDist + 3;
      numchunks  = 0;
    }
    numchunks += 1;
    param.chunk_dist = hamDist;
    param.numblocks  = numBlocks;
    param.numchunks  = numchunks;
  } while (calcMissingEdgeRatio(param) >= _missingratio);
}

void SketchSort::zNormalization() {
  size_t num = fvs.size();
  for (size_t d = 0; d < dim_; ++d) {
    double ave = 0.0;
    for (size_t i = 0; i < num; ++i) 
      ave += fvs[i][d]; 
    ave /= (double)num;

    double sigma = 0.0;
    for (size_t i = 0; i < num; ++i) 
      sigma += (fvs[i][d] - ave) * (fvs[i][d] - ave);
    sigma /= (num - 1);
    sigma = sqrt(sigma);
    if (sigma < 1.0e-30)
      continue;
    for (size_t i = 0; i < num; ++i) 
      fvs[i][d] = (fvs[i][d] - ave)/sigma;
  }
}

void SketchSort::MinMaxNormalization() {
  size_t num = fvs.size();
  for (size_t d = 0; d < dim_; ++d) {
    double minVal = DBL_MAX, maxVal = -DBL_MAX;
    for (size_t i = 0; i < num; ++i) {
      if (fvs[i][d] > maxVal)
	maxVal = fvs[i][d];
      if (fvs[i][d] < minVal)
	minVal = fvs[i][d];
    }
    if (minVal < 1.0e-30 && maxVal < 1.0e-30)
      continue;
    for (size_t i = 0; i < num; ++i) 
      fvs[i][d] = (fvs[i][d] - minVal)/(maxVal - minVal);
  }
}

void SketchSort::run(const char *fname, const char *oname, 
		     uint32_t _numblocks, 
		     uint32_t _dist, 
		     float    _minmaxDist,
		     uint32_t _numchunks, 
		     bool     _autoFlag,
		     float    _missingratio,
		     bool     _zNormalization, 
		     bool     _minmaxNormalization) 
{
  params param;
  param.numblocks  = _numblocks;
  param.numchunks  = _numchunks;
  param.chunk_dist = _dist;
  param.minmaxDist = _minmaxDist;
  param.projectDim = 32;

  if (_autoFlag) {
    cerr << "turn parameters such that the missing edge ratio is no more than " << _missingratio << endl;
    turnParameters(_missingratio, param);
    cout << "turned parameters:" << endl;
    cout << "hamming distance threshold: " << param.chunk_dist << endl;
    cout << "number of blocks: " << param.numblocks << endl;
    cout << "number of chunks: "  << param.numchunks << endl;
    cout << endl;
  }

  std::ofstream ofs(oname);
  param.os = &ofs;

  cout << "missing edge ratio:" << calcMissingEdgeRatio(param) << endl;

  cerr << "start reading" << endl;
  double readstart = clock();
  readFeature(fname);
  double readend   = clock();
  cerr << "end reading" << endl;
  cout << "readtime:" << (readend - readstart)/(double)CLOCKS_PER_SEC << endl;
  
  num_char = 256;
  double totalstart = clock();

  if (_zNormalization)  {
    cout << "z-normalizing" << endl;
    zNormalization();
  }
  else if (_minmaxNormalization) {
    cout << "minmax normalizing" << endl;
    MinMaxNormalization();
  }

  param.counter = new uint32_t[num_char];
  
  cout << "number of data:" << fvs.size() << endl;
  cout << "length of strings per chunk:" << param.projectDim << endl;
  cout << "number of chunks:" << param.numchunks << endl;
  cout << "total length of strings:" << param.projectDim * param.numchunks << endl;

  double projectstart = clock();
  cerr << "start projection" << endl;
  vector<uint8_t*> sig;
  projectVectors(param.projectDim * param.numchunks, sig, param);
  //read(fname, sig, param);
  cerr << "end projection" << std::endl;
  double projectend = clock();
  cout << "projecttime:" << (projectend - projectstart)/(double)CLOCKS_PER_SEC << endl;
  param.pchunks = new pstat[param.numchunks + 1];
  for (int i = 1; i <= (int)param.numchunks; i++) {
    param.pchunks[i].start = (int)ceil((double)param.seq_len*((double)(i - 1)/(double)param.numchunks)) + 1;
    param.pchunks[i].end   = (int)ceil((double)param.seq_len*(double)i/(double)param.numchunks);
  }

  double msmtime = 0.0;

  cerr << "chunk distance:" << param.chunk_dist << endl;
  cerr << "the number of blocks:" << param.numblocks << endl;
  param.pos = new pstat[param.numblocks + 1];
  for (int i = 1; i <= (int) param.numchunks; i++) {
    param.chunk_len   = param.pchunks[i].end - param.pchunks[i].start + 1;
    param.start_chunk = param.pchunks[i].start; 
    param.end_chunk   = param.pchunks[i].end;
    param.cchunk      = i;
    for (int j = 1; j <= (int)param.numblocks; j++) {
      param.pos[j].start = (int)ceil((double)param.chunk_len*((double)(j - 1)/(double)param.numblocks)) + param.pchunks[i].start;
      param.pos[j].end   = (int)ceil((double)param.chunk_len*(double)j/(double)param.numblocks) + param.pchunks[i].start - 1;
    }
    cerr << "start enumeration chunk no " << i << endl;
    double msmstart = clock();
    multi_classification(sig, 1, 0, param.num_seq - 1, param);
    double msmend   = clock();
    msmtime += (msmend - msmstart)/(double)CLOCKS_PER_SEC;
  }
  cout << "msmtime:" << msmtime << endl;

  double totalend = clock();
  cout << "cputime:" << (totalend - totalstart)/(double)CLOCKS_PER_SEC << endl;

  ofs.close();
  // destructor
  delete p;
  delete[] param.counter;
  delete[] param.pchunks;
  delete[] param.pos;

  return;
}
