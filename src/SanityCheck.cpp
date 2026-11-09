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
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"

namespace {
// why a merged count finds no usable list, see Scores::countMergedTargets
const char* const kUnusableListReasons =
    "the scores of a test bin cannot be rescaled, or targets and decoys "
    "separate too well in a test bin (test, --static), in the merged list or "
    "in its best PSM per spectrum (-Y) or per peptide";
}  // namespace

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
                                      test_fdr, reportedLists_));
  }
  // The second chance counts the vectors of a reset only when it runs;
  // without --init-weights they are the initial ones, and training changes
  // w, so they are kept for it.
  initWeights_.clear();
  if (useSecondChance(w.size())) {
    initWeights_ = w;
  }
  return initPositives_;
}

void SanityCheck::getDefaultDirection(const vector<Scores>& trainset,
                                      vector<vector<double> >& w) const {
    
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

bool SanityCheck::wrongSignOfDefaultDirection(
    const vector<vector<double> >& w) const {
  if (!initDefaultDir_) {
    return false;
  }
  for (size_t set = 0; set < w.size(); ++set) {
    if (w[set][static_cast<std::size_t>(abs(initDefaultDir_) - 1)] *
            initDefaultDir_ <= 0) {
      return true;
    }
  }
  return false;
}

vector<vector<double> > SanityCheck::resetVectors(
    const vector<Scores>& trainset) const {
  if (initWeightFN_.empty()) {
    return initWeights_;
  }
  vector<vector<double> > w(
      trainset.size(), vector<double>(FeatureNames::getNumFeatures() + 1));
  // the reset itself reports the features it chooses
  const int verbosity = VERB;
  Globals::getInstance()->setVerbose(std::min(verbosity, 1));
  try {
    getDefaultDirection(trainset, w);
  } catch (const MyException&) {
    w.clear();
  } catch (...) {
    Globals::getInstance()->setVerbose(verbosity);
    throw;
  }
  Globals::getInstance()->setVerbose(verbosity);
  return w;
}

bool SanityCheck::secondChance(const vector<Scores>& testset,
                               const vector<Scores>& trainset,
                               const vector<vector<double> >& w,
                               int& trainedCount) const {
  // the anchors of the training sets exist for test bins of any size
  trainedCount = Scores::countMergedTargets(
      testset, trainset, w, MergeNormalization::TRAIN_DECOYS, mergeFdr_,
      skipNormalizeScores_, test_fdr_, reportedLists_);
  // The fallback that the second chance replaces, counted as the run would
  // report it: each fold's default vector, merged with the run's scheme
  // (test up to three folds) and counted at the test FDR. A fallback whose
  // vectors cannot be chosen, or whose list cannot be used, would end the
  // run and counts 0.
  int resetCount = 0;
  const vector<vector<double> > resetWeights = resetVectors(trainset);
  if (resetWeights.size() == w.size()) {
    resetCount = std::max(
        0, Scores::countMergedTargets(testset, trainset, resetWeights,
                                      mergeScheme_, mergeFdr_,
                                      skipNormalizeScores_, test_fdr_,
                                      reportedLists_));
  }
  const bool countsKeep = trainedCount > 0 && trainedCount >= resetCount;
  // A run that keeps the trained SVMs is merged with the scheme of
  // setSecondChanceMergeScheme. Another scheme than train-decoys, given by
  // --merge-normalization, can end that run where train-decoys does not:
  // test checks the separation of each bin, and the scores of test, train
  // and none are merged on another scale. Its list therefore has to be
  // usable as well; the count of train-decoys still decides. --static
  // merges the raw scores, which the count above has already checked.
  const bool otherKeptScheme =
      !skipNormalizeScores_ &&
      secondChanceMergeScheme_ != MergeNormalization::TRAIN_DECOYS;
  const bool unusableAsKept =
      countsKeep && otherKeptScheme &&
      Scores::countMergedTargets(testset, trainset, w,
                                 secondChanceMergeScheme_, mergeFdr_,
                                 skipNormalizeScores_, test_fdr_,
                                 reportedLists_) < 0;
  // the sign check of --default-direction applies as after a passed count
  const bool wrongSign =
      countsKeep && !unusableAsKept && wrongSignOfDefaultDirection(w);
  const bool keep = countsKeep && !unusableAsKept && !wrongSign;
  cerr << "Second chance on the merged out-of-fold list"
       << (skipNormalizeScores_ ? "" : ", rescaled by the training decoys")
       << ": ";
  if (trainedCount < 0) {
    cerr << "the trained SVMs give no usable list (" << kUnusableListReasons
         << ")";
  } else {
    cerr << trainedCount << " targets with q<" << test_fdr_
         << " for the trained SVMs, " << resetCount
         << " for the default vectors of a reset";
    if (!skipNormalizeScores_) {
      cerr << " merged with " << mergeNormalizationName(mergeScheme_);
    }
    if (unusableAsKept) {
      cerr << "; merged with "
           << mergeNormalizationName(secondChanceMergeScheme_)
           << ", as a run that keeps them would be, the trained SVMs give no "
           << "usable list (" << kUnusableListReasons << ")";
    }
  }
  if (overRule_) {
    cerr << "; --override keeps the trained SVMs in any case." << endl;
  } else if (keep) {
    cerr << "; the trained SVMs are kept";
    if (!skipNormalizeScores_) {
      cerr << ", merged with "
           << mergeNormalizationName(secondChanceMergeScheme_);
    }
    cerr << "." << endl;
  } else if (wrongSign) {
    cerr << "; the trained SVMs are not kept, as the weight of the main "
         << "scoring direction (--default-direction) has the wrong sign."
         << endl;
  } else {
    cerr << "; the trained SVMs are not kept." << endl;
  }
  return keep;
}

bool SanityCheck::validateDirection(vector<Scores>& testset,
                                    vector<Scores>& trainset,
                                    vector<vector<double> >& w) {
  secondChanceKept_ = false;
  const bool pooled = usePooledCount(w.size());
  int overFDR = 0;
  if (pooled) {
    // on copies: the test bins are only scored once the direction is decided
    overFDR = Scores::countMergedTargets(testset, trainset, w, mergeScheme_,
                                         mergeFdr_, skipNormalizeScores_,
                                         test_fdr_, reportedLists_);
    if (overFDR < 0) {
      cerr << "The trained weights give no usable merged out-of-fold list: "
           << kUnusableListReasons << endl;
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
  bool keptBySecondChance = false;
  if (overFDR <= 0 || initPositives_ > overFDR) {
    if (overFDR <= 0) {
      cerr << "No targets found with q<" << test_fdr_ << endl;
    } else {
      cerr << "Less identifications (" << overFDR << " vs " << initPositives_
          << ") after percolator processing than before processing" << endl;
    }
    // Where the per-bin count would reset, the second chance judges the
    // trained SVMs on the merged out-of-fold list, which small test bins
    // do not limit. It decides nothing under --override.
    const bool secondChanceTried = useSecondChance(w.size());
    int pooledCount = 0;
    keptBySecondChance = secondChanceTried &&
                         secondChance(testset, trainset, w, pooledCount) &&
                         !overRule_;
    if (!keptBySecondChance) {
      resetDirection(trainset, w);
      if (overFDR <= 0 && !overRule_) {
        // not for an unusable merged list (-1), whose reason is printed
        if (pooled || (secondChanceTried && pooledCount == 0)) {
          cerr << "Hint: q<" << test_fdr_ << " needs about "
               << 1.0 / test_fdr_ << " targets above the best-scoring decoy "
               << "of the merged list; a larger --testFDR is an alternative."
               << endl;
        } else if (w.size() > 1u && !secondChanceTried) {
          cerr << "Hint: each of the " << w.size() << " test folds may be too "
               << "small to reach q<" << test_fdr_ << ". --num-folds 1 trains "
               << "and evaluates on the same PSMs," << endl
               << "so its q values and PEPs are optimistic, and wrong for "
               << "inputs without signal; a larger --testFDR is an "
               << "alternative." << endl;
        }
      }
      return false;
    }
  }
  if (wrongSignOfDefaultDirection(w)) {
    if (VERB > 1) {
      cerr << "Warning, wrong sign of the weight for main scoring direction"
           << endl;
    }
    resetDirection(trainset, w);
    return false;
  }
  if (pooled) {
    // the direction is kept: score the test bins with their folds' weights
    // for the merge, as the count of the bins check does
    for (size_t set = 0; set < w.size(); ++set) {
      testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  // the per-bin count has already scored the test bins with w
  secondChanceKept_ = keptBySecondChance;
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

