/*
 * Copyright 2021 Brian Raiter <breadbox@muppetlabs.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <map>
#include "SetHandler.h"
#include "CrossValidation.h"
#include "Globals.h"

/* A subclass of CrossValidation that gives us access to some
 * protected fields.
 */
class CrossValidationEx : public CrossValidation {
public:
    CrossValidationEx(bool quickValidation,
                      bool reportPerformanceEachIteration, double testFdr,
                      double selectionFdr, double initialSelectionFdr,
                      double selectedCpos, double selectedCneg,
                      unsigned int niter, bool usePi0,
                      unsigned int nestedXvalBins, bool trainBestPositive,
                      unsigned int numThreads, bool skipNormalizeScores, 
                      double decoyFractionTraining, unsigned int numFolds)
        : CrossValidation(quickValidation, reportPerformanceEachIteration,
                          testFdr, selectionFdr, initialSelectionFdr,
                          selectedCpos, selectedCneg, niter, usePi0,
                          nestedXvalBins, trainBestPositive, numThreads,
                          skipNormalizeScores, decoyFractionTraining, numFolds) { }
    virtual ~CrossValidationEx() { }
    int doStepEx(Normalizer* pNorm, double selectionFdr) {
        return doStep(pNorm, selectionFdr);
    }
    std::vector< std::vector<double> > const& weights(void) const {
        return weights_;
    }
    std::vector<Scores> const& testBins(void) const {
        return testScores_;
    }
};

/*
 * Unit tests for the CrossValidation class.
 *
 * Note: These tests take over 250 ms to complete (on my computer).
 * They are okay here as long as these are the only slow tests in the
 * unit test suite. But at some point we will probably want to move
 * this class into a separate suite of integration tests.
 */

class CrossValidationTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    virtual void TearDown();
    int origVerbose;
};

void CrossValidationTest::SetUp()
{
    origVerbose = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(0);
}

void CrossValidationTest::TearDown()
{
    Globals::getInstance()->setVerbose(origVerbose);
}

// With numFolds folds each PSM is in numFolds - 1 training sets, and the
// estimated number of positives is their sum over the training sets divided
// by numFolds - 1; with one fold the only training set holds all PSMs. The
// estimate is therefore about the number of positives in the whole data set
// (N) for any number of folds.
TEST_F(CrossValidationTest, doStepTest)
{
  // the default 3 folds first, then other numbers of folds; with 2 folds the
  // training sets hold only about N/2 targets above the first decoy, too few
  // to reach q < testFdr in some of them, and the estimate is about N/2
  for (unsigned int numFolds : {3u, 1u, 5u, 10u}) {
    SCOPED_TRACE(numFolds);
    // Note that we use an elevated testFdr (0.02 instead of 0.01),
    // to compensate for the small size of the data sets.
    int const N = 100;
    double const testFdr = 0.02;

    CrossValidationEx *crossValidation =
            new CrossValidationEx(false,   // quickValidation
                                  false,   // reportEachIteration
                                  testFdr, // testFdr
                                  0.01,    // selectionFdr
                                  0.01,    // initialSectionFdr
                                  0.0,     // selectedCpos
                                  0.0,     // selectedCneg
                                  10,      // nIter
                                  true,    // usePi0
                                  1,       // nestedXvalBins
                                  false,   // trainBestPositive
                                  1,       // numThreads
                                  false,   // skipNormalizeScores
                                  1.0,     // decoyFractionTraining
                                  numFolds);

    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    SetHandler setHandler(0);
    int scanNumber = 1;

    // Our data set has two features. One flips between 0 and 1
    // without regard for label. The other is consistently increasing,
    // but with the targets' values slightly less than the decoys'
    // values.

    PsmStore psmStore;
    DataSet *targets = new DataSet();
    DataSet *decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    for (int i = 0 ; i < 2 * N ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = new double[2];
        psm->features[0] = static_cast<double>(i) / N - 1E-6;
        psm->features[1] = static_cast<double>(i % 2);
        psm->scan = scanNumber++;
        targets->registerPsm(psm);
    }
    for (int i = 0 ; i < N ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = new double[2];
        psm->features[0] = static_cast<double>(i) / N;
        psm->features[1] = static_cast<double>(i % 2);
        psm->scan = scanNumber++;
        decoys->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(targets));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(decoys));

    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    setHandler.normalizeFeatures(*norm);
    Normalizer* pNorm = norm.get();

    Scores scores(true);
    scores.populateWithPSMs(setHandler);

    // The pre-iteration setup should select the first feature to
    // be weighted.

    SanityCheck check;
    int numInit = crossValidation->preIterationSetup(scores, check, pNorm);
    std::vector< std::vector<double> > const& w = crossValidation->weights();
    ASSERT_EQ(numFolds, w.size());
    for (unsigned int i = 0 ; i < numFolds ; ++i) {
        EXPECT_EQ(1.0, w[i][0]);
        EXPECT_EQ(0.0, w[i][1]);
    }

    // One step of the training algorithm should find at least N positives.
    // The number can be a bit higher because some targets with score < N are
    // counted as positives. This can be either because the highest scoring decoy 
    // in the cross validation fold has a lower score, or because it narrowly passes 
    // the testFDR threshold.
    // The estimate is truncated, and the trained direction can move single
    // targets below a decoy in some training sets, so it can be slightly
    // below N for other numbers of folds.
    int estimatedNumPositives = crossValidation->doStepEx(pNorm, 0.01);
    if (numFolds == 3u) {
      EXPECT_LE(N, estimatedNumPositives);
    } else {
      EXPECT_LE(N * (1.0 - testFdr), estimatedNumPositives);
    }
    EXPECT_GE(N * (1.0 + testFdr), estimatedNumPositives);

    delete crossValidation;
  }
  FeatureNames::resetNumFeatures();
}

// The direction check of trainAndMerge
struct DirectionCheckOptions {
    bool overrule = false;  // --override
    double testFdr = 0.05;  // a tiny one makes the check fail
};

// Trains and merges with numFolds folds and the given merge scheme on 450
// targets and 450 decoys. Half of the targets score higher in both features,
// so that the SVM improves on the initial direction of a single feature.
// Returns the merged score of each scan and the final weights, and checks
// that the final weights of each fold give the merged scores of its test
// bin: each PSM was scored by the fold of its own test bin.
static void trainAndMerge(unsigned int numFolds, MergeNormalization scheme,
                          std::map<unsigned int, double>& scoreOf,
                          std::vector<std::vector<double> >& weights,
                          std::vector<double>& features,
                          const DirectionCheckOptions& checkOptions =
                              DirectionCheckOptions()) {
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    const int numPsms = 900;
    features.assign(2 * numPsms, 0.0);
    SetHandler setHandler(0);
    PsmStore psmStore;
    DataSet* targets = new DataSet();
    DataSet* decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    unsigned int seed = 2024u;
    for (int i = 0; i < numPsms; ++i) {
        seed = seed * 1103515245u + 12345u;
        const bool isTarget = i % 2 == 0;
        PSMDescription* psm = psmStore.add(PSMDescription());
        psm->features = &features[2 * i];
        const bool isCorrect = isTarget && i % 4 == 0;
        psm->features[0] = static_cast<double>((seed >> 8) % 1000) / 250.0 +
                           (isCorrect ? 2.5 : 0.0);
        psm->features[1] = static_cast<double>((seed >> 4) % 1000) / 250.0 +
                           (isCorrect ? 2.5 : 0.0);
        psm->scan = static_cast<unsigned int>(1 + i);
        (isTarget ? targets : decoys)->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(targets));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(decoys));
    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    setHandler.normalizeFeatures(*norm);

    PseudoRandom::setSeed(1u);
    CrossValidationEx crossValidation(false, false, checkOptions.testFdr,
                                      0.05, 0.05, 0.0, 0.0, 10, false, 1,
                                      false, 1, false, 1.0, numFolds);
    crossValidation.setMergeNormalization(scheme);
    Scores scores(false);
    scores.populateWithPSMs(setHandler);
    SanityCheck check;
    check.setOverrule(checkOptions.overrule);
    crossValidation.preIterationSetup(scores, check, norm.get());
    crossValidation.train(norm.get());
    crossValidation.postIterationProcessing(scores, check);

    scoreOf.clear();
    for (const ScoreHolder& sh : scores) {
        scoreOf[sh.pPSM->scan] = sh.score;
    }
    weights = crossValidation.weights();
    std::map<const PSMDescription*, std::size_t> binOf;
    for (std::size_t bin = 0; bin < numFolds; ++bin) {
        for (const ScoreHolder& sh : crossValidation.testBins()[bin]) {
            binOf[sh.pPSM] = bin;
        }
    }
    ASSERT_EQ(static_cast<std::size_t>(numPsms), binOf.size());
    ASSERT_EQ(static_cast<unsigned int>(numPsms), scores.size());
    for (const ScoreHolder& sh : scores) {
        const double expected =
            scores.calcScore(sh.pPSM->features, weights[binOf[sh.pPSM]]);
        EXPECT_NEAR(expected, sh.score,
                    1e-12 * std::max(1.0, std::fabs(sh.score)))
            << "bin " << binOf[sh.pPSM];
    }
    FeatureNames::resetNumFeatures();
}

// trainAndMerge finds each test bin scored with its own fold's final
// weights. When the check passes, they are the trained ones. When it fails,
// here by a test FDR that no target reaches, they are, with --override, the
// trained ones, which were not trained on that bin, and otherwise the
// default vector of a single feature chosen on the fold's training set.
// Before, all bins were scored with the first fold's vector after a failed
// check; after a reset this failed even where the folds chose the same
// feature, since each fold's default vector kept its trained bias.
TEST_F(CrossValidationTest, FailedDirectionCheckScoresEachBinWithItsOwnVector)
{
    auto isSingleFeature = [](const std::vector<double>& w) {
        return (w[0] == 0.0) != (w[1] == 0.0);
    };
    std::map<unsigned int, double> scoreOf;
    std::vector<double> features;
    for (unsigned int numFolds : {3u, 5u}) {
        {
            SCOPED_TRACE(::testing::Message() << numFolds
                                              << " folds, check passes");
            std::vector<std::vector<double> > w;
            trainAndMerge(numFolds, MergeNormalization::TEST, scoreOf, w,
                          features);
            ASSERT_EQ(numFolds, w.size());
            for (const std::vector<double>& wFold : w) {
                EXPECT_FALSE(isSingleFeature(wFold));  // trained, not reset
            }
        }
        for (bool overrule : {false, true}) {
            SCOPED_TRACE(::testing::Message() << numFolds
                                              << " folds, failed check, "
                                              << "override " << overrule);
            DirectionCheckOptions failing;
            failing.overrule = overrule;
            failing.testFdr = 1e-6;
            std::vector<std::vector<double> > w;
            trainAndMerge(numFolds, MergeNormalization::TEST, scoreOf, w,
                          features, failing);
            ASSERT_EQ(numFolds, w.size());
            for (const std::vector<double>& wFold : w) {
                EXPECT_EQ(!overrule, isSingleFeature(wFold));
            }
            if (overrule) {
                // each fold keeps its own trained weights
                for (unsigned int fold = 1; fold < numFolds; ++fold) {
                    EXPECT_NE(w[0], w[fold]);
                }
            }
        }
    }
}

// With one fold the training set is the test bin, so the training-set
// anchors of --merge-normalization train are the test-bin anchors: the
// merged scores and the weights equal those of the default scheme exactly.
// With three folds every scheme gives scores that the final weights
// describe; none leaves the SVM weights as trained.
TEST_F(CrossValidationTest, MergeNormalizationSchemes)
{
    std::map<unsigned int, double> scoreTest, scoreOther;
    std::vector<std::vector<double> > wTest, wOther;
    std::vector<double> featuresTest, featuresOther;
    {
        SCOPED_TRACE("one fold");
        trainAndMerge(1u, MergeNormalization::TEST, scoreTest, wTest,
                      featuresTest);
        trainAndMerge(1u, MergeNormalization::TRAIN, scoreOther, wOther,
                      featuresOther);
        EXPECT_EQ(wTest, wOther);
        ASSERT_EQ(scoreTest.size(), scoreOther.size());
        ASSERT_EQ(900u, scoreTest.size());
        // one PSM per scan, compared in scan order
        auto it = scoreTest.begin();
        auto jt = scoreOther.begin();
        for (; it != scoreTest.end(); ++it, ++jt) {
            EXPECT_EQ(it->second, jt->second);
        }
    }
    for (MergeNormalization scheme :
         {MergeNormalization::TEST, MergeNormalization::TRAIN,
          MergeNormalization::TRAIN_DECOYS, MergeNormalization::NONE}) {
        SCOPED_TRACE(::testing::Message() << "three folds, scheme "
                                          << static_cast<int>(scheme));
        trainAndMerge(3u, scheme, scoreOther, wOther, featuresOther);
        EXPECT_EQ(900u, scoreOther.size());
        EXPECT_EQ(3u, wOther.size());
    }
}

// More folds than spectra leave a fold without PSMs, which
// preIterationSetup reports as an error. With as many folds as spectra of
// equal size, each fold gets one spectrum.
TEST_F(CrossValidationTest, EmptyFoldIsAnError)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    SetHandler setHandler(0);
    PsmStore psmStore;
    DataSet *targets = new DataSet();
    DataSet *decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    std::vector<double> features(2 * 8);
    // 2 spectra with 2 targets and 2 decoys each
    for (int i = 0 ; i < 8 ; ++i) {
        PSMDescription *psm = psmStore.add(PSMDescription());
        psm->features = &features[2 * i];
        psm->features[0] = static_cast<double>(i);
        psm->features[1] = static_cast<double>(i % 3);
        psm->scan = 1 + i / 4;
        (i % 2 ? decoys : targets)->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(targets));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(decoys));
    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::STDV);
    setHandler.normalizeFeatures(*norm);

    for (unsigned int numFolds : {1u, 2u, 3u}) {
        SCOPED_TRACE(numFolds);
        CrossValidation crossValidation(false, false, 0.01, 0.01, 0.01, 0.0,
                                        0.0, 10, false, 1, false, 1, false,
                                        1.0, numFolds);
        Scores scores(false);
        scores.populateWithPSMs(setHandler);
        EXPECT_EQ(2u, scores.getNumSpectra());
        SanityCheck check;
        if (numFolds <= 2u) {
            EXPECT_NO_THROW(
                crossValidation.preIterationSetup(scores, check, norm.get()));
        } else {
            EXPECT_THROW(
                crossValidation.preIterationSetup(scores, check, norm.get()),
                MyException);
        }
    }
    FeatureNames::resetNumFeatures();
}
