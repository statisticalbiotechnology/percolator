/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"

void SanityCheck::checkAndSetDefaultDir() {
  if (!initDefaultDir_ && initDefaultDirName_.size() > 0) {
    int sign = 1;
    if (initDefaultDirName_[0] == '-') {
      initDefaultDirName_.erase(0,1);
      sign = -1;
    }
    initDefaultDir_ = sign * DataSet::getFeatureNames().getFeatureNumber(initDefaultDirName_);
    if (initDefaultDir_ == 0) {
      ostringstream temp;
      temp << "ERROR: Initial direction feature name \"" << initDefaultDirName_ << "\" not found" << std::endl;
      throw MyException(temp.str());
    }
  }
}

int SanityCheck::getInitDirection(vector<Scores>& testset,
                                  vector<Scores>& trainset,
                                  const Normalizer* pNorm,
                                  vector<vector<double> > & w,
                                  double test_fdr,
                                  double initial_train_fdr) {
  test_fdr_= test_fdr;
  initial_train_fdr_ = initial_train_fdr;
  if (initWeightFN_.size() > 0) {
    vector<double> ww(FeatureNames::getNumFeatures() + 1);
    ifstream weightStream(initWeightFN_.data(), ios::in);
    if (weightStream.is_open()) {
      readWeights(weightStream, ww);
      weightStream.close();
      assert(pNorm);
      pNorm->normalizeweight(ww, w[0]);
      for (size_t set = 1; set < w.size(); ++set) {
        copy(w[0].begin(), w[0].end(), w[set].begin());
      }
    } else {
      std::cerr << "WARNING: Could not find weights input file " << initWeightFN_
                << ". Using default weights instead." << std::endl;
      getDefaultDirection(trainset, w);
    }
  } else {
    getDefaultDirection(trainset, w);
  }
  initPositives_ = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    initPositives_ += testset[set].calcScoresAndQvals(w[set], test_fdr);
  }
  // The pooled check compares the trained direction with the initial one
  // counted in the same way: per fold, merged and pooled. The bins are
  // still scored above, so that both checks leave them in the same state.
  if (usePooledCount(w.size())) {
    initPositives_ = std::max(
        0, Scores::countMergedTargets(testset, trainset, w, mergeScheme_,
                                      mergeFdr_, skipNormalizeScores_,
                                      test_fdr));
  }
  return initPositives_;
}

void SanityCheck::getDefaultDirection(vector<Scores>& trainset,
                                      vector<vector<double> >& w) {
    
  //If I have not been given a initial direction
  if (!initDefaultDir_) {
    //if the default weights from the input file are not present
    if (defaultWeights_.size() == 0) {
      // Set init direction to be the most discriminative direction
      for (size_t set = 0; set < w.size(); ++set) {
        if (VERB > 1) std::cerr << "Split " << set + 1 << ":\t";
        trainset[set].getInitDirection(initial_train_fdr_, w[set]);
      }
    } else {
      // I want to assign the default vector that is present in the input file
      for (size_t set = 0; set < w.size(); ++set) {
        for (size_t ix = 0; ix < w[set].size(); ix++) {
          w[set][ix] = 0;
          if (ix < defaultWeights_.size()){
              w[set][ix] = defaultWeights_[ix];
          }          
        }
      }
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      for (size_t ix = 0; ix < w[set].size(); ix++) {
        w[set][ix] = 0;
      }
      w[set][static_cast<std::size_t>(abs(initDefaultDir_) - 1)] = (initDefaultDir_ < 0 ? -1 : 1);
    }
  }
}

bool SanityCheck::validateDirection(vector<Scores>& testset,
                                    vector<Scores>& trainset,
                                    vector<vector<double> >& w) {
  const bool pooled = usePooledCount(w.size());
  int overFDR = 0;
  if (pooled) {
    // on copies: the test bins are only scored once the direction is decided
    overFDR = Scores::countMergedTargets(testset, trainset, w, mergeScheme_,
                                         mergeFdr_, skipNormalizeScores_,
                                         test_fdr_);
    if (overFDR < 0) {
      cerr << "The scores of a test bin cannot be rescaled to merge the "
           << "bins with the trained weights" << endl;
      resetDirection(trainset, w);
      return false;
    }
    if (VERB > 1) {
      cerr << "Found " << overFDR << " targets with q<" << test_fdr_
           << " in the merged out-of-fold list, " << initPositives_
           << " in the initial direction" << endl;
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      overFDR += testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  if (overFDR <= 0) {
    cerr << "No targets found with q<" << test_fdr_ << endl;
    resetDirection(trainset, w);
    if (pooled && !overRule_) {
      cerr << "Hint: q<" << test_fdr_ << " needs about " << 1.0 / test_fdr_
           << " targets above the best-scoring decoy of the merged list; "
           << "a larger --testFDR is an alternative." << endl;
    } else if (w.size() > 1u && !overRule_) {
      cerr << "Hint: each of the " << w.size() << " test folds may be too "
           << "small to reach q<" << test_fdr_ << ". --num-folds 1 trains and "
           << "evaluates on the same PSMs," << endl
           << "so its q values and PEPs are optimistic, and wrong for inputs "
           << "without signal; a larger --testFDR is an alternative." << endl;
    }
    return false;
  }
  if (initPositives_ > overFDR) {
    cerr << "Less identifications (" << overFDR << " vs " << initPositives_
        << ") after percolator processing than before processing" << endl;
    resetDirection(trainset, w);
    return false;
  }
  if (initDefaultDir_) {
    for (size_t set = 0; set < w.size(); ++set) {
      if (w[set][static_cast<std::size_t>(abs(initDefaultDir_) - 1)] * initDefaultDir_ <= 0) {
        if (VERB > 1) {
          cerr
            << "Warning, wrong sign of the weight for main scoring direction"
            << endl;
        }
        resetDirection(trainset, w);
        return false;
      }
    }
  }
  if (pooled) {
    // the direction is kept: score the test bins with their folds' weights
    // for the merge, as the count of the bins check does
    for (size_t set = 0; set < w.size(); ++set) {
      testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  return true;
}

void SanityCheck::readWeights(istream& weightStream, vector<double>& w) {
  std::string line;
  // Find first line with weights, but skip this and get second line containing raw features
  while (std::getline(weightStream, line)) {
    if (line[0] == '-' || (line[0] >= '0' && line[0] <= '9')) break;
  }
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream >> w[ix];
  }
  cerr << "Read weights from file" << endl;
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    cerr << w[ix] << "\t";
  }
  cerr << endl;
}

void SanityCheck::resetDirection(vector<Scores>& trainset,
                                 vector<vector<double> >& w) {
  if (!overRule_) {
    cerr << "Resetting score vector, using default vector. Use --override flag to prevent this." << endl;
    getDefaultDirection(trainset, w);
  }
}

