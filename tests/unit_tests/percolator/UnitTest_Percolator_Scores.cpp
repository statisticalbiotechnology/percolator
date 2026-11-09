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

/*
 * Unit tests for the Scores and ScoreHolder classes.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <utility>
#include "SetHandler.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "PosteriorEstimator.h"
#include "SanityCheck.h"
#include "Scores.h"

// Some strings in alphabetical order.
static std::string const psmNames[] = { "ABC", "DEF", "GHI", "JKL", "MNO" };

class ScoreHolderTest : public ::testing::Test {
  protected:
    bool checkOrder(std::vector<ScoreHolder> const *scores, ...);
};

// A simple function that verifies that the scores in a vector are in
// a specific order. (The function is variadic so that the caller
// doesn't need to declare an extra array.)
bool ScoreHolderTest::checkOrder(std::vector<ScoreHolder> const *scores, ...)
{
    va_list values;
    va_start(values, scores);
    for (std::vector<ScoreHolder>::const_iterator it = scores->begin() ;
             it != scores->end() ; ++it) {
        double score = va_arg(values, double);
        if (it->score != score)
            return false;
    }
    va_end(values);
    return true;
}

// Verify the ScoreHolder's ordering functions.
TEST_F(ScoreHolderTest, CheckOrderingFunctions)
{
    std::vector<ScoreHolder> scores;

    // peptides are assigned in alphabetic order
    // scan values are assigned [ 2, 3, 4, 0, 1 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
        pPSM->scan = (i + 2) % 5;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    std::sort(scores.begin(), scores.end(), lexicOrderProb());
    ASSERT_TRUE(checkOrder(&scores, 1.0, 2.0, 3.0, 4.0, 5.0));
    std::sort(scores.begin(), scores.end(), OrderScanMassCharge());
    ASSERT_TRUE(checkOrder(&scores, 4.0, 5.0, 1.0, 2.0, 3.0));
    
    // specFileNr values are assigned [ 3, 3, 4, 4, 5]
    // scan values are assigned [ 0, 0, 1, 1, 2]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
        pPSM->specFileNr = 3 + i / 2;
        pPSM->scan = i / 2;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanHash());
    for (int i = 0 ; i < 5 ; ++i) {
        std::cerr << scores.at(i).score << std::endl;
    }
    ASSERT_TRUE(checkOrder(&scores, 3.0, 4.0, 1.0, 2.0, 5.0));

    // specFileNr values are assigned [ 0, 0, 1, 1, 2]
    // scan values are assigned [ 2, 3, 4, 0, 1 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
        pPSM->specFileNr = i / 2;
        pPSM->scan = (i + 2) % 5;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanMassCharge());
    ASSERT_TRUE(checkOrder(&scores, 1.0, 2.0, 4.0, 3.0, 5.0));

    // scan values are assigned [ 2, 2, 1, 1, 0 ]
    // labels are assigned [ -1, +1, -1, +1, -1 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription();
        pPSM->scan = 2 - i / 2;
        scores.push_back(ScoreHolder(1.0 + i, (i % 2 ? LabelType::TARGET : LabelType::DECOY), pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanLabel());
    ASSERT_TRUE(checkOrder(&scores, 5.0, 4.0, 3.0, 2.0, 1.0));

    // scan values are assigned [ 2, 2, 1, 1, 0 ]
    // labels are assigned [ -1, +1, -1, +1, -1 ]
    // specFileNr values are assigned [ 0, 0, 1, 1, 2]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription();
        pPSM->specFileNr = i / 2;
        pPSM->scan = 2 - i / 2;
        scores.push_back(ScoreHolder(1.0 + i, (i % 2 ? LabelType::TARGET : LabelType::DECOY), pPSM));
    }
    std::sort(scores.begin(), scores.end(), OrderScanLabel());
    ASSERT_TRUE(checkOrder(&scores, 2.0, 1.0, 4.0, 3.0, 5.0));
}

// Verify the ScoreHolder's uniqueness functions.
TEST_F(ScoreHolderTest, CheckUniquenessFilter)
{
    std::vector<ScoreHolder> scores;

    // use unique peptides but repeating scan values
    // peptides are assigned in alphabetic order
    // scan values are assigned [ 0, 0, 1, 1, 2 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
        pPSM->scan = i / 2;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    ASSERT_EQ(5, scores.size());
    scores.erase(std::unique(scores.begin(), scores.end(), UniqueScanLabel()),
                 scores.end());
    ASSERT_EQ(3, scores.size());

    // use unique peptides but repeating specFileNr and scan values
    // peptides are assigned in alphabetic order
    // specFileNr values are assigned [ 0, 0, 0, 1, 1]
    // scan values are assigned [ 0, 0, 1, 1, 2 ]
    scores.clear();
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *pPSM = new PSMDescription(psmNames[i]);
        pPSM->specFileNr = i / 3;
        pPSM->scan = i / 2;
        scores.push_back(ScoreHolder(1.0 + i, LabelType::TARGET, pPSM));
    }
    ASSERT_EQ(5, scores.size());
    scores.erase(std::unique(scores.begin(), scores.end(), UniqueScanLabel()),
                 scores.end());
    ASSERT_EQ(4, scores.size());
}


class ScoresTest : public ::testing::Test {
  protected:
    virtual void SetUp();
    virtual void TearDown();
  private:
    int origVerbose;
};

void ScoresTest::SetUp()
{
    origVerbose = Globals::getInstance()->getVerbose();
    Globals::getInstance()->setVerbose(0);
}

void ScoresTest::TearDown()
{
    Globals::getInstance()->setVerbose(origVerbose);
}

// Test the basic populateWithPsms() method.
TEST_F(ScoresTest, CheckPopulating)
{
    Scores scores(true);
    SetHandler setHandler(0);
    PsmStore psmStore;
    DataSet *set1 = new DataSet();
    DataSet *set2 = new DataSet();
    set1->setLabel(LabelType::TARGET);
    set2->setLabel(LabelType::DECOY);
    for (int i = 0 ; i < 5 ; ++i) {
        PSMDescription *psm;
        psm = psmStore.add(PSMDescription(psmNames[i]));
        psm->scan = 100 + i;
        set1->registerPsm(psm);
        psm = psmStore.add(PSMDescription(psmNames[i]));
        psm->scan = 105 + i;
        set2->registerPsm(psm);
    }
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set1));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set2));

    scores.populateWithPSMs(setHandler);
    ASSERT_EQ(10, scores.size());

    unsigned scanValue = 100;
    for (std::vector<ScoreHolder>::const_iterator it = scores.begin() ;
            it != scores.end() ;
            it++, scanValue++) {
        ASSERT_EQ(scanValue, it->pPSM->scan);
    }
}

// Test that populateWithPsms() refuses empty sets.
TEST_F(ScoresTest, CheckPopulatingEmpty)
{
    Scores scores(true);
    SetHandler setHandler(0);
    EXPECT_THROW(scores.populateWithPSMs(setHandler), MyException);
    DataSet *set1 = new DataSet();
    DataSet *set2 = new DataSet();
    set1->setLabel(LabelType::TARGET);
    set2->setLabel(LabelType::DECOY);
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set1));
    setHandler.push_back_dataset(std::unique_ptr<DataSet>(set2));
    EXPECT_THROW(scores.populateWithPSMs(setHandler), MyException);
}

// createXvalSetsBySpectrum with N folds puts every PSM in the test set of
// exactly one fold and in the training sets of the other N - 1 folds, keeps
// the PSMs of a spectrum in one fold and gives folds of about equal size.
// With one fold, the training set and the test set both hold all PSMs.
TEST_F(ScoresTest, CreateXvalSetsBySpectrumWithNFolds)
{
    // 2 files with 37 spectra each, with 1 to 4 PSMs per spectrum and the
    // same scan numbers in both files
    std::vector<PSMDescription> psms;
    psms.reserve(400);
    for (unsigned int file = 0; file < 2; ++file) {
        for (unsigned int scan = 0; scan < 37; ++scan) {
            for (unsigned int k = 0; k <= (scan * 7 + file) % 4; ++k) {
                psms.emplace_back();
                psms.back().specFileNr = file;
                psms.back().scan = 1000 + scan;
            }
        }
    }
    const std::size_t maxPsmsPerSpectrum = 4;
    ASSERT_LT(psms.size(), 400u);  // no reallocation
    auto spectrumOf = [](const PSMDescription* psm) {
        return std::make_pair(psm->specFileNr, psm->scan);
    };

    for (unsigned int numFolds : {1u, 2u, 3u, 5u, 10u}) {
        Scores scores(false);
        for (std::size_t i = 0; i < psms.size(); ++i) {
            scores.addScoreHolder(ScoreHolder(
                0.0, i % 3 ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
        }
        scores.recalculateSizes();
        ASSERT_EQ(74u, scores.getNumSpectra());

        std::vector<Scores> train, test;
        scores.createXvalSetsBySpectrum(train, test, numFolds);
        ASSERT_EQ(numFolds, train.size());
        ASSERT_EQ(numFolds, test.size());

        std::map<const PSMDescription*, unsigned int> testFold, numTrain;
        std::map<std::pair<unsigned int, unsigned int>, unsigned int>
            spectrumFold;
        std::size_t numTest = 0;
        for (unsigned int fold = 0; fold < numFolds; ++fold) {
            EXPECT_EQ(test[fold].size(),
                      static_cast<unsigned int>(
                          std::distance(test[fold].begin(), test[fold].end())));
            EXPECT_LE(test[fold].size(),
                      (psms.size() + numFolds - 1) / numFolds +
                          maxPsmsPerSpectrum - 1)
                << numFolds << " folds";
            EXPECT_GT(test[fold].size(), 0u) << numFolds << " folds";
            for (const ScoreHolder& sh : test[fold]) {
                ++numTest;
                EXPECT_EQ(0u, testFold.count(sh.pPSM));  // in one test set
                testFold[sh.pPSM] = fold;
                auto spectrum = spectrumOf(sh.pPSM);
                if (spectrumFold.count(spectrum)) {
                    EXPECT_EQ(spectrumFold[spectrum], fold);
                } else {
                    spectrumFold[spectrum] = fold;
                }
            }
        }
        EXPECT_EQ(psms.size(), numTest);
        EXPECT_EQ(74u, spectrumFold.size());
        for (unsigned int fold = 0; fold < numFolds; ++fold) {
            for (const ScoreHolder& sh : train[fold]) {
                ++numTrain[sh.pPSM];
                // a PSM is not trained on in the fold it is tested in, unless
                // there is only one fold
                if (numFolds > 1) {
                    EXPECT_NE(testFold[sh.pPSM], fold);
                }
            }
        }
        EXPECT_EQ(psms.size(), numTrain.size());
        for (const auto& psmCount : numTrain) {
            EXPECT_EQ(std::max(1u, numFolds - 1), psmCount.second);
        }
        if (numFolds == 1) {
            EXPECT_EQ(psms.size(), train[0].size());
            EXPECT_EQ(psms.size(), test[0].size());
            EXPECT_EQ(scores.posSize(), train[0].posSize());
            EXPECT_EQ(scores.negSize(), test[0].negSize());
        }
    }
}

// The test fold of each PSM as createXvalSetsBySpectrum assigned it before
// the fold capacity was fixed, from the same random draws: remain counts
// PSMs, but the folds are filled with whole spectra. sorted holds the PSMs in
// the order of the split (OrderScanHash); labels receives the labels of the
// test sets, pseudo targets included.
static std::vector<std::size_t> previousFoldAssignment(
    const Scores& sorted, unsigned int numFolds, double decoyFractionTraining,
    std::vector<LabelType>& labels) {
    std::vector<int> remain(numFolds);
    int fold = static_cast<int>(numFolds);
    int ix = static_cast<int>(std::distance(sorted.begin(), sorted.end()));
    while (fold--) {
        remain[static_cast<std::size_t>(fold)] = ix / (fold + 1);
        ix -= remain[static_cast<std::size_t>(fold)];
    }
    std::vector<std::size_t> foldOf;
    labels.clear();
    unsigned int previousScan = sorted.begin()->pPSM->scan;
    std::size_t randIndex = PseudoRandom::lcg_rand() % numFolds;
    for (const ScoreHolder& sh : sorted) {
        if (sh.pPSM->scan != previousScan) {
            randIndex = PseudoRandom::lcg_rand() % numFolds;
            while (remain[randIndex] <= 0) {
                randIndex = PseudoRandom::lcg_rand() % numFolds;
            }
        }
        LabelType label = sh.label;
        if (numFolds > 1u && decoyFractionTraining < 1.0 &&
            label == LabelType::DECOY &&
            PseudoRandom::lcg_uniform_rand() > decoyFractionTraining) {
            label = LabelType::PSEUDO_TARGET;
        }
        foldOf.push_back(randIndex);
        labels.push_back(label);
        --remain[randIndex];
        previousScan = sh.pPSM->scan;
    }
    return foldOf;
}

// The test fold and label of each PSM after createXvalSetsBySpectrum, in the
// order it leaves scores in, and the next random number.
struct FoldSplit {
    std::vector<std::size_t> foldOf;
    std::vector<LabelType> labels;
    unsigned long nextRandom = 0;
    std::size_t numEmptyFolds = 0;
};

static FoldSplit splitIntoFolds(Scores& scores, unsigned int numFolds,
                                double decoyFractionTraining,
                                unsigned long seed) {
    PseudoRandom::setSeed(seed);
    std::vector<Scores> train, test;
    scores.createXvalSetsBySpectrum(train, test, numFolds,
                                    decoyFractionTraining);
    FoldSplit split;
    split.nextRandom = PseudoRandom::lcg_rand();
    std::map<const PSMDescription*, std::pair<std::size_t, LabelType> > of;
    for (std::size_t fold = 0; fold < numFolds; ++fold) {
        if (test[fold].begin() == test[fold].end()) {
            ++split.numEmptyFolds;
        }
        for (const ScoreHolder& sh : test[fold]) {
            of[sh.pPSM] = std::make_pair(fold, sh.label);
        }
    }
    for (const ScoreHolder& sh : scores) {
        split.foldOf.push_back(of.at(sh.pPSM).first);
        split.labels.push_back(of.at(sh.pPSM).second);
    }
    return split;
}

// The same with the assignment before the fix, on the order that
// createXvalSetsBySpectrum left scores in.
static FoldSplit previousSplit(const Scores& sorted, unsigned int numFolds,
                               double decoyFractionTraining,
                               unsigned long seed) {
    PseudoRandom::setSeed(seed);
    FoldSplit split;
    split.foldOf = previousFoldAssignment(sorted, numFolds,
                                          decoyFractionTraining, split.labels);
    split.nextRandom = PseudoRandom::lcg_rand();
    std::vector<bool> used(numFolds, false);
    split.numEmptyFolds = numFolds;
    for (std::size_t fold : split.foldOf) {
        if (!used[fold]) {
            used[fold] = true;
            --split.numEmptyFolds;
        }
    }
    return split;
}

// The first index where two splits differ in the fold or the label of a
// PSM, or their size if they agree. A plain loop: under GCC's parallel mode,
// which the build uses, comparing vectors with == starts threads.
static std::size_t firstDifference(const FoldSplit& a, const FoldSplit& b) {
    const std::size_t n = std::min(a.foldOf.size(), b.foldOf.size());
    for (std::size_t i = 0; i < n; ++i) {
        if (a.foldOf[i] != b.foldOf[i] || a.labels[i] != b.labels[i]) {
            return i;
        }
    }
    return n;
}

// The spectra as the split delimits them: a new one starts where the scan
// number changes between neighbours.
static std::size_t numSplitSpectra(const Scores& sorted) {
    std::size_t numSpectra = 0;
    const PSMDescription* previous = NULL;
    for (const ScoreHolder& sh : sorted) {
        if (previous == NULL || previous->scan != sh.pPSM->scan) {
            ++numSpectra;
        }
        previous = sh.pPSM;
    }
    return numSpectra;
}

// PSMs of a separate search: each spectrum alternates targets and decoys.
static void fillScores(Scores& scores, std::vector<PSMDescription>& psms) {
    for (std::size_t i = 0; i < psms.size(); ++i) {
        psms[i].setLineNr(static_cast<unsigned int>(i + 1));
        scores.addScoreHolder(ScoreHolder(
            0.0, i % 2 ? LabelType::DECOY : LabelType::TARGET, &psms[i]));
    }
    scores.recalculateSizes();
}

// The folds are filled with whole spectra, but their capacity is counted in
// PSMs. With uneven spectra and many folds, a fold could take more spectra
// than its share and leave another fold without any, which ended the run
// (the empty-fold error of CrossValidation::preIterationSetup). With at
// least as many spectra as folds, every fold now gets one, with and without
// the pseudo targets of RESET; the assignment before the fix left folds
// empty on these inputs.
TEST_F(ScoresTest, CreateXvalSetsBySpectrumFillsEveryFold)
{
    // 120 spectra of 1 to 5 PSMs
    std::vector<PSMDescription> psms;
    for (unsigned int scan = 1; scan <= 120; ++scan) {
        for (unsigned int k = 0; k <= (scan * 7) % 5; ++k) {
            psms.emplace_back();
            psms.back().scan = scan;
        }
    }
    std::size_t previousEmpty = 0;
    for (unsigned int numFolds : {120u, 119u, 100u, 60u, 30u, 10u}) {
        for (double decoyFractionTraining : {1.0, 0.5}) {
            SCOPED_TRACE(::testing::Message() << numFolds << " folds, "
                                              << decoyFractionTraining);
            Scores scores(false);
            fillScores(scores, psms);
            const FoldSplit split =
                splitIntoFolds(scores, numFolds, decoyFractionTraining, 1u);
            EXPECT_EQ(120u, numSplitSpectra(scores));
            EXPECT_EQ(0u, split.numEmptyFolds);
            // a spectrum is in one fold
            for (std::size_t i = 1; i < psms.size(); ++i) {
                if (scores.begin()[i].pPSM->scan ==
                    scores.begin()[i - 1].pPSM->scan) {
                    EXPECT_EQ(split.foldOf[i - 1], split.foldOf[i]);
                }
            }
            previousEmpty += previousSplit(scores, numFolds,
                                           decoyFractionTraining, 1u)
                                 .numEmptyFolds;
        }
    }
    EXPECT_LT(0u, previousEmpty);
}

// The fix changes no assignment that left every fold with a spectrum, nor
// the random draws: on random inputs (1 to 6 PSMs per spectrum, one or two
// spectrum files with shared scan numbers, up to 3 folds more than spectra,
// with and without RESET's pseudo targets) the folds, labels and the next
// random number equal those of the assignment before the fix whenever that
// left no fold empty, and also when there are fewer spectra than folds.
// Otherwise no fold is empty now.
TEST_F(ScoresTest, CreateXvalSetsBySpectrumKeepsAssignmentsWithoutEmptyFolds)
{
    std::mt19937 generator(404u);
    std::size_t numUnchanged = 0, numFewerSpectra = 0, numFilled = 0;
    for (int run = 0; run < 2000; ++run) {
        const unsigned int numScans = 1u + generator() % 60u;
        const unsigned int maxPsms = 1u + generator() % 6u;
        const unsigned int numFiles = 1u + generator() % 2u;
        std::vector<PSMDescription> psms;
        for (unsigned int file = 0; file < numFiles; ++file) {
            for (unsigned int scan = 1; scan <= numScans; ++scan) {
                if (numFiles > 1u && generator() % 3u == 0u) {
                    continue;  // not every scan in every file
                }
                const unsigned int n = 1u + generator() % maxPsms;
                for (unsigned int k = 0; k < n; ++k) {
                    psms.emplace_back();
                    psms.back().specFileNr = file;
                    psms.back().scan = 1000u + scan;
                }
            }
        }
        if (psms.empty()) {
            continue;
        }
        Scores scores(false);
        fillScores(scores, psms);
        const unsigned int numSpectra =
            static_cast<unsigned int>(scores.getNumSpectra());
        const unsigned int numFolds =
            1u + generator() % std::min(numSpectra + 3u, 100u);
        const double decoyFractionTraining =
            generator() % 2u ? 1.0 : 0.5;
        const unsigned long seed = 1u + generator() % 1000u;
        SCOPED_TRACE(::testing::Message()
                     << "run " << run << ": " << psms.size() << " PSMs, "
                     << numSpectra << " spectra, " << numFolds << " folds, "
                     << decoyFractionTraining << ", seed " << seed);
        const FoldSplit split =
            splitIntoFolds(scores, numFolds, decoyFractionTraining, seed);
        const FoldSplit previous =
            previousSplit(scores, numFolds, decoyFractionTraining, seed);
        if (previous.numEmptyFolds == 0u ||
            numSplitSpectra(scores) < numFolds) {
            ASSERT_EQ(previous.foldOf.size(), split.foldOf.size());
            EXPECT_EQ(previous.foldOf.size(), firstDifference(previous, split));
            EXPECT_EQ(previous.nextRandom, split.nextRandom);
            EXPECT_EQ(previous.numEmptyFolds, split.numEmptyFolds);
            ++(previous.numEmptyFolds == 0u ? numUnchanged : numFewerSpectra);
        } else {
            EXPECT_EQ(0u, split.numEmptyFolds);
            ++numFilled;
        }
    }
    EXPECT_LT(500u, numUnchanged);
    EXPECT_LT(50u, numFewerSpectra);
    EXPECT_LT(50u, numFilled);
}

// Scores::getInitDirection picks the single feature that separates best and
// sets every other weight to 0, the bias included. After a failed direction
// check, SanityCheck::resetDirection calls it on each fold's training set
// with the fold's trained weights in place, which kept their bias.
TEST_F(ScoresTest, DefaultDirectionAfterResetHasZeroBias)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    // two folds of 20 targets and 20 decoys; feature 1 separates the
    // targets (higher values) from the decoys, feature 2 does not
    const std::size_t numFolds = 2, numPsms = 40;
    std::vector<double> features(2 * numFolds * numPsms);
    std::vector<PSMDescription> psms(numFolds * numPsms);
    std::vector<Scores> testBins(numFolds, Scores(false)),
        trainBins(numFolds, Scores(false));
    for (std::size_t i = 0; i < psms.size(); ++i) {
        const bool isTarget = i % 2 == 0;
        psms[i].features = &features[2 * i];
        psms[i].features[0] =
            static_cast<double>(i % numPsms) + (isTarget ? 100.0 : 0.0);
        psms[i].features[1] = static_cast<double>((i * 7) % 11);
        psms[i].scan = static_cast<unsigned int>(i + 1);
        const ScoreHolder sh(0.0,
                             isTarget ? LabelType::TARGET : LabelType::DECOY,
                             &psms[i]);
        const std::size_t fold = i / numPsms;
        testBins[fold].addScoreHolder(sh);
        trainBins[1 - fold].addScoreHolder(sh);
    }
    for (std::size_t fold = 0; fold < numFolds; ++fold) {
        testBins[fold].recalculateSizes();
        trainBins[fold].recalculateSizes();
    }
    const std::vector<double> defaultVector = {1.0, 0.0, 0.0};
    {
        std::vector<double> direction = {0.3, -2.0, 7.5};
        trainBins[0].getInitDirection(0.01, direction);
        EXPECT_EQ(defaultVector, direction);
    }
    for (bool overrule : {false, true}) {
        SCOPED_TRACE(overrule);
        SanityCheck check;
        check.setDirectionCheck(SanityCheck::DirectionCheck::BINS);
        check.setOverrule(overrule);
        std::vector<std::vector<double> > w(numFolds,
                                            std::vector<double>(3, 0.0));
        // a test FDR that no target reaches makes the check fail
        check.getInitDirection(testBins, trainBins, NULL, w, 1e-6, 0.01);
        const std::vector<std::vector<double> > trained = {{0.5, 0.25, 3.0},
                                                           {0.4, -0.1, -2.0}};
        w = trained;
        EXPECT_FALSE(check.validateDirection(testBins, trainBins, w));
        if (overrule) {
            EXPECT_EQ(trained, w);
        } else {
            for (const std::vector<double>& wFold : w) {
                EXPECT_EQ(defaultVector, wFold);
            }
        }
    }
    FeatureNames::resetNumFeatures();
}

// normalizeScores rescales the scores of each cross-validation bin to a
// common scale before the bins are merged. When the median decoy does not
// score below the score at the FDR threshold, which happens for tiny inputs,
// this is impossible and an error, unless the bin is the only one: its raw
// SVM scores and weights are then kept. A bin that can be rescaled is
// rescaled in the same way with and without singleBin.
TEST_F(ScoresTest, NormalizeScoresKeepsRawScoresOfASingleBin)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    const double fdr = 0.01;
    const std::vector<double> rawWeights = {0.5, -0.25, 1.0};
    std::vector<PSMDescription> psms(20);
    for (std::size_t i = 0; i < psms.size(); ++i) {
        psms[i].scan = static_cast<unsigned int>(i);
    }
    // a decoy above three targets, sorted by score; with a single decoy
    // there is no median decoy below the top score
    const double unscalableScores[] = {4.0, 3.0, 2.0, 1.0};
    auto unscalable = [&psms, &unscalableScores]() {
        Scores scores(false);
        for (std::size_t i = 0; i < 4; ++i) {
            scores.addScoreHolder(ScoreHolder(unscalableScores[i],
                i == 0 ? LabelType::DECOY : LabelType::TARGET, &psms[i]));
        }
        scores.recalculateSizes();
        return scores;
    };
    auto expectUnchanged = [&unscalableScores](const Scores& scores) {
        ASSERT_EQ(4u, scores.size());
        std::size_t i = 0;
        for (const ScoreHolder& sh : scores) {
            EXPECT_EQ(unscalableScores[i++], sh.score);
        }
    };

    {
        SCOPED_TRACE("normalizeScores");
        Scores scores = unscalable();
        scores.calcQvals(fdr);
        std::vector<double> w = rawWeights;
        EXPECT_THROW(scores.normalizeScores(fdr, w), MyException);
        EXPECT_NO_THROW(scores.normalizeScores(fdr, w, true));
        EXPECT_EQ(rawWeights, w);
        expectUnchanged(scores);
    }
    {
        SCOPED_TRACE("merge of one bin");
        std::vector<Scores> bins(1, unscalable());
        std::vector<std::vector<double> > w(1, rawWeights);
        Scores merged(false);
        EXPECT_NO_THROW(merged.merge(bins, fdr, false, w));
        EXPECT_EQ(rawWeights, w[0]);
        expectUnchanged(merged);
    }
    {
        SCOPED_TRACE("merge of two bins");
        std::vector<Scores> bins(2, unscalable());
        std::vector<std::vector<double> > w(2, rawWeights);
        Scores merged(false);
        EXPECT_THROW(merged.merge(bins, fdr, false, w), MyException);
    }
    {
        SCOPED_TRACE("rescalable bin");
        // 10 targets above 10 decoys
        Scores scores[2] = {Scores(false), Scores(false)};
        std::vector<double> w[2] = {rawWeights, rawWeights};
        for (int k = 0; k < 2; ++k) {
            for (std::size_t i = 0; i < 20; ++i) {
                scores[k].addScoreHolder(ScoreHolder(
                    20.0 - static_cast<double>(i),
                    i < 10 ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
            }
            scores[k].recalculateSizes();
            scores[k].calcQvals(fdr);
            EXPECT_NO_THROW(scores[k].normalizeScores(fdr, w[k], k == 1));
        }
        EXPECT_NE(rawWeights, w[0]);
        EXPECT_EQ(w[0], w[1]);
        ASSERT_EQ(20u, scores[0].size());
        auto it0 = scores[0].begin();
        auto it1 = scores[1].begin();
        EXPECT_NE(20.0, it0->score);
        for (; it0 != scores[0].end(); ++it0, ++it1) {
            EXPECT_EQ(it0->score, it1->score);
        }
    }
    FeatureNames::resetNumFeatures();
}

// Synthetic PSMs with two features for the merge-normalization tests: a
// third of the targets score like the decoys, the others higher.
class MergeNormalizationData {
  public:
    static const std::size_t numFeatures = 2;
    explicit MergeNormalizationData(int numPsms)
        : features(numPsms * numFeatures), psms(numPsms), labels(numPsms) {
        unsigned int seed = 4711u;
        for (int i = 0; i < numPsms; ++i) {
            seed = seed * 1103515245u + 12345u;
            bool isTarget = (seed >> 16) % 2 == 0;
            double* row = &features[i * numFeatures];
            row[0] = static_cast<double>((seed >> 8) % 101) / 10.0 +
                     (isTarget && i % 3 ? 4.0 : 0.0);
            row[1] = static_cast<double>((seed >> 4) % 13) / 4.0;
            psms[i].features = row;
            psms[i].scan = static_cast<unsigned int>(i);
            labels[i] = isTarget ? LabelType::TARGET : LabelType::DECOY;
        }
    }
    // bin b holds the PSMs i with i % numBins == b, scored with w[b]
    std::vector<Scores> bins(std::size_t numBins, bool usePi0,
                             std::vector<std::vector<double> > w) {
        std::vector<Scores> result(numBins, Scores(usePi0));
        for (std::size_t i = 0; i < psms.size(); ++i) {
            result[i % numBins].addScoreHolder(
                ScoreHolder(0.0, labels[i], &psms[i]));
        }
        for (std::size_t b = 0; b < numBins; ++b) {
            result[b].recalculateSizes();
            result[b].calcScores(w[b]);
        }
        return result;
    }
    std::size_t binOf(const ScoreHolder& sh, std::size_t numBins) const {
        return static_cast<std::size_t>(sh.pPSM - &psms[0]) % numBins;
    }
    std::vector<double> features;
    std::vector<PSMDescription> psms;
    std::vector<LabelType> labels;
};

// With the TRAIN scheme each bin is rescaled by the anchors of its fold's
// training set. When the training set is the test bin itself, as with one
// fold, the merge gives exactly the scores and weights of the TEST scheme,
// also with the pi0 of separate searches and with several bins. With every
// scheme the rescaled weights still give the merged scores, and with
// TRAIN_DECOYS the decoys of a bin that is its own training set get mean 0
// and standard deviation 1.
TEST_F(ScoresTest, MergeWithTrainingAnchorsEqualsTestWhenTrainIsTest)
{
    typedef MergeNormalizationData Data;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(Data::numFeatures);
    const double fdr = 0.05;
    Data data(600);
    const std::vector<std::vector<double> > binWeights = {
        {1.0, -0.5, 0.25}, {2.0, 0.3, -1.0}, {0.7, -0.1, 0.0}};
    const MergeNormalization schemes[] = {
        MergeNormalization::TEST, MergeNormalization::TRAIN,
        MergeNormalization::TRAIN_DECOYS, MergeNormalization::NONE};

    for (bool usePi0 : {false, true}) {
        for (std::size_t numBins : {1u, 3u}) {
            SCOPED_TRACE(::testing::Message() << "usePi0 " << usePi0
                                              << ", " << numBins << " bins");
            const std::vector<std::vector<double> > w0(
                binWeights.begin(), binWeights.begin() + numBins);
            std::map<MergeNormalization, std::vector<std::vector<double> > >
                weights;
            std::map<MergeNormalization, Scores> merged;
            for (MergeNormalization scheme : schemes) {
                std::vector<Scores> bins = data.bins(numBins, usePi0, w0);
                // the training set of each bin is the bin itself
                std::vector<Scores> trainSets = data.bins(numBins, usePi0, w0);
                std::vector<ScoreScale> scales;
                std::vector<std::vector<double> > w = w0;
                for (std::size_t b = 0; b < numBins; ++b) {
                    scales.push_back(
                        trainSets[b].calcScoreScale(w[b], scheme, fdr));
                }
                Scores m(usePi0);
                ASSERT_NO_THROW(m.merge(bins, fdr, false, w, scheme,
                                        scheme == MergeNormalization::TEST
                                            ? std::vector<ScoreScale>()
                                            : scales));
                // the weights give the merged scores
                for (const ScoreHolder& sh : m) {
                    const double expected =
                        m.calcScore(sh.pPSM->features,
                                    w[data.binOf(sh, numBins)]);
                    EXPECT_NEAR(expected, sh.score,
                                1e-12 * std::max(1.0, std::fabs(expected)));
                }
                if (scheme == MergeNormalization::NONE) {
                    EXPECT_EQ(w0, w);
                } else {
                    EXPECT_NE(w0, w);
                }
                weights[scheme] = w;
                merged.insert(std::make_pair(scheme, m));
            }

            const Scores& test = merged.at(MergeNormalization::TEST);
            const Scores& train = merged.at(MergeNormalization::TRAIN);
            EXPECT_EQ(weights[MergeNormalization::TEST],
                      weights[MergeNormalization::TRAIN]);
            ASSERT_EQ(test.size(), train.size());
            for (auto it = test.begin(), jt = train.begin(); it != test.end();
                 ++it, ++jt) {
                EXPECT_EQ(it->pPSM, jt->pPSM);
                EXPECT_EQ(it->score, jt->score);
            }

            if (numBins == 1u) {
                double sum = 0.0, sumSquares = 0.0;
                std::size_t numDecoys = 0u;
                for (const ScoreHolder& sh :
                     merged.at(MergeNormalization::TRAIN_DECOYS)) {
                    if (sh.isDecoy()) {
                        sum += sh.score;
                        sumSquares += sh.score * sh.score;
                        ++numDecoys;
                    }
                }
                ASSERT_GT(numDecoys, 100u);
                const double mean = sum / numDecoys;
                EXPECT_NEAR(0.0, mean, 1e-12);
                EXPECT_NEAR(1.0, (sumSquares - numDecoys * mean * mean) /
                                     (numDecoys - 1), 1e-12);
            }
        }
    }
    FeatureNames::resetNumFeatures();
}

// The NONE scheme merges the raw SVM scores and keeps the weights, also for
// bins that are too small to be rescaled. Training-set anchors that do not
// exist are an error with several bins, as test-bin anchors are, while a
// single bin keeps its raw scores.
TEST_F(ScoresTest, MergeWithoutNormalizationKeepsRawScores)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    const double fdr = 0.01;
    const std::vector<double> rawWeights = {0.5, -0.25, 1.0};
    std::vector<PSMDescription> psms(8);
    // two bins without a median decoy below the score at q < fdr
    const double rawScores[] = {4.0, 3.0, 2.0, 1.5, 1.0, 10.0, 0.5, 0.2};
    const bool isDecoy[] = {true, false, false, true, false,
                            false, true, false};
    const std::size_t binEnd[] = {5, 8};
    auto makeBins = [&](std::size_t numBins) {
        std::vector<Scores> bins(numBins, Scores(false));
        for (std::size_t i = 0; i < binEnd[numBins - 1]; ++i) {
            bins[i < binEnd[0] ? 0 : 1].addScoreHolder(ScoreHolder(
                rawScores[i], isDecoy[i] ? LabelType::DECOY : LabelType::TARGET,
                &psms[i]));
        }
        for (Scores& bin : bins) {
            bin.recalculateSizes();
        }
        return bins;
    };
    std::map<const PSMDescription*, double> rawScoreOf;
    for (std::size_t i = 0; i < psms.size(); ++i) {
        rawScoreOf[&psms[i]] = rawScores[i];
    }

    {
        SCOPED_TRACE("none");
        std::vector<Scores> bins = makeBins(2);
        std::vector<std::vector<double> > w(2, rawWeights);
        Scores merged(false);
        ASSERT_NO_THROW(
            merged.merge(bins, fdr, false, w, MergeNormalization::NONE));
        EXPECT_EQ(std::vector<std::vector<double> >(2, rawWeights), w);
        ASSERT_EQ(psms.size(), merged.size());
        double previous = 1e9;
        for (const ScoreHolder& sh : merged) {
            EXPECT_EQ(rawScoreOf[sh.pPSM], sh.score);
            EXPECT_LE(sh.score, previous);  // sorted
            previous = sh.score;
        }
    }
    {
        SCOPED_TRACE("test");
        std::vector<Scores> bins = makeBins(2);
        std::vector<std::vector<double> > w(2, rawWeights);
        Scores merged(false);
        EXPECT_THROW(merged.merge(bins, fdr, false, w), MyException);
    }
    for (MergeNormalization scheme :
         {MergeNormalization::TRAIN, MergeNormalization::TRAIN_DECOYS}) {
        SCOPED_TRACE(static_cast<int>(scheme));
        const ScoreScale missing;  // scale 0: the anchors do not exist
        {
            std::vector<Scores> bins = makeBins(2);
            std::vector<std::vector<double> > w(2, rawWeights);
            Scores merged(false);
            EXPECT_THROW(merged.merge(bins, fdr, false, w, scheme,
                                      std::vector<ScoreScale>(2, missing)),
                         MyException);
            // the scales of all bins are needed
            EXPECT_THROW(merged.merge(bins, fdr, false, w, scheme,
                                      std::vector<ScoreScale>(1)),
                         MyException);
        }
        {
            std::vector<Scores> bins = makeBins(1);
            std::vector<std::vector<double> > w(1, rawWeights);
            Scores merged(false);
            EXPECT_NO_THROW(merged.merge(bins, fdr, false, w, scheme,
                                         std::vector<ScoreScale>(1, missing)));
            EXPECT_EQ(rawWeights, w[0]);
            ASSERT_EQ(binEnd[0], merged.size());
            for (const ScoreHolder& sh : merged) {
                EXPECT_EQ(rawScoreOf[sh.pPSM], sh.score);
            }
        }
    }
    FeatureNames::resetNumFeatures();
}

// A hand-made example of the pooled direction check: three test bins of six
// PSMs with a single feature x, scored by three folds whose weights put x on
// different scales: x, 2x + 1 and x / 2 - 1. The decoys of every fold's
// training set have x = -1, 0 and 1, so train-decoys maps each fold's
// scores back to x exactly, and the merged out-of-fold list is ordered by x:
//   x    10 9.8 9.5 9 8.8 8.5 8 7.8 7.5 7 6.8 6.5 6.2 6 5.5 5 3 1
//   T/D  T  T   T   T T   T   T T   T   T T   T   D   D D   T T D
// With q = (decoys + 1) / targets the first 12 targets have q = 1/12 and the
// targets at 5 and 3 have q = 4/14. Within each bin only 4 targets precede
// the best decoy, which gives q = 1/4 at best.
class PooledCountExample {
  public:
    static constexpr std::size_t numBins = 3, binSize = 6, trainSize = 5;
    PooledCountExample()
        : features(numBins * (binSize + trainSize)),
          psms(features.size()),
          labels(features.size()) {
        const double testX[numBins][binSize] = {{10, 9, 8, 7, 6, 5},
                                                {9.5, 8.5, 7.5, 6.5, 5.5, 1},
                                                {9.8, 8.8, 7.8, 6.8, 6.2, 3}};
        const bool testIsTarget[numBins][binSize] = {{1, 1, 1, 1, 0, 1},
                                                     {1, 1, 1, 1, 0, 0},
                                                     {1, 1, 1, 1, 0, 1}};
        const double trainX[trainSize] = {-1, 0, 1, 4, 5};
        const bool trainIsTarget[trainSize] = {0, 0, 0, 1, 1};
        std::size_t i = 0;
        for (std::size_t bin = 0; bin < numBins; ++bin) {
            for (std::size_t j = 0; j < binSize; ++j, ++i) {
                set(i, testX[bin][j], testIsTarget[bin][j]);
            }
            for (std::size_t j = 0; j < trainSize; ++j, ++i) {
                set(i, trainX[j], trainIsTarget[j]);
            }
        }
    }
    // bin b, or the training set of fold b, unscored
    std::vector<Scores> testBins(std::size_t n = numBins) {
        return sets(n, 0, binSize);
    }
    std::vector<Scores> trainBins(std::size_t n = numBins) {
        return sets(n, binSize, trainSize);
    }
    // the score of each test PSM rescaled by hand with the mean and the
    // sample standard deviation of its fold's training decoys, or raw
    std::vector<std::pair<double, bool> > handMerged(bool rescale) const {
        std::vector<std::pair<double, bool> > merged;
        for (std::size_t bin = 0; bin < numBins; ++bin) {
            const std::size_t first = bin * (binSize + trainSize);
            double sum = 0.0, sumSquares = 0.0;
            int numDecoys = 0;
            for (std::size_t i = first + binSize; i < first + binSize + trainSize;
                 ++i) {
                if (labels[i] == LabelType::DECOY) {
                    sum += score(bin, i);
                    ++numDecoys;
                }
            }
            const double mean = sum / numDecoys;
            for (std::size_t i = first + binSize; i < first + binSize + trainSize;
                 ++i) {
                if (labels[i] == LabelType::DECOY) {
                    sumSquares += (score(bin, i) - mean) * (score(bin, i) - mean);
                }
            }
            const double sd = std::sqrt(sumSquares / (numDecoys - 1));
            for (std::size_t i = first; i < first + binSize; ++i) {
                merged.push_back(std::make_pair(
                    rescale ? (score(bin, i) - mean) / sd : score(bin, i),
                    labels[i] == LabelType::TARGET));
            }
        }
        return merged;
    }
    std::vector<std::vector<double> > w = {
        {1.0, 0.0}, {2.0, 1.0}, {0.5, -1.0}};
    std::vector<double> features;
    std::vector<PSMDescription> psms;
    std::vector<LabelType> labels;

  private:
    void set(std::size_t i, double x, bool isTarget) {
        features[i] = x;
        psms[i].features = &features[i];
        psms[i].scan = static_cast<unsigned int>(i);
        labels[i] = isTarget ? LabelType::TARGET : LabelType::DECOY;
    }
    double score(std::size_t bin, std::size_t i) const {
        return w[bin][0] * features[i] + w[bin][1];
    }
    std::vector<Scores> sets(std::size_t n, std::size_t offset,
                             std::size_t size) {
        std::vector<Scores> result(n, Scores(false));
        for (std::size_t bin = 0; bin < n; ++bin) {
            const std::size_t first = bin * (binSize + trainSize) + offset;
            for (std::size_t i = first; i < first + size; ++i) {
                result[bin].addScoreHolder(ScoreHolder(0.0, labels[i], &psms[i]));
            }
            result[bin].recalculateSizes();
        }
        return result;
    }
};

// The targets with q = (decoys + 1) / targets < testFdr, counted directly.
static int directCount(std::vector<std::pair<double, bool> > merged,
                       double testFdr) {
    std::sort(merged.begin(), merged.end(),
              [](const std::pair<double, bool>& a,
                 const std::pair<double, bool>& b) {
                  return a.first > b.first;
              });
    std::vector<double> q(merged.size());
    int decoys = 0, targets = 0;
    for (std::size_t i = 0; i < merged.size(); ++i) {
        (merged[i].second ? targets : decoys) += 1;
        q[i] = (decoys + 1.0) / std::max(1, targets);
    }
    for (std::size_t i = merged.size() - 1; i-- > 0;) {
        q[i] = std::min(q[i], q[i + 1]);
    }
    int count = 0;
    for (std::size_t i = 0; i < merged.size(); ++i) {
        if (merged[i].second && q[i] < testFdr) {
            ++count;
        }
    }
    return count;
}

// The pooled count of Scores::countMergedTargets equals a direct count of
// the q values of the hand-made merged list, and the count of merge and
// calcQvals, which give the final list; it leaves the bins as they are. A
// training set whose decoys do not spread cannot rescale its bin (-1), and
// with one bin the pooled count is the bin's own count.
TEST_F(ScoresTest, CountMergedTargetsEqualsDirectCount)
{
    typedef PooledCountExample Example;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(1);
    Example example;
    const double fdr = 0.01;  // the anchors of train-decoys do not use it
    const MergeNormalization trainDecoys = MergeNormalization::TRAIN_DECOYS;

    // the hand-made count of the comment above
    EXPECT_EQ(12, directCount(example.handMerged(true), 0.25));
    EXPECT_EQ(14, directCount(example.handMerged(true), 0.3));
    for (double testFdr : {0.25, 0.3, 0.5}) {
        SCOPED_TRACE(testFdr);
        std::vector<Scores> bins = example.testBins();
        const std::vector<Scores> trainBins = example.trainBins();
        const int count = Scores::countMergedTargets(
            bins, trainBins, example.w, trainDecoys, fdr, false, testFdr);
        EXPECT_EQ(directCount(example.handMerged(true), testFdr), count);
        // on copies: the bins are not scored
        for (const Scores& bin : bins) {
            for (const ScoreHolder& sh : bin) {
                EXPECT_EQ(0.0, sh.score);
            }
        }
        // the final list of the same bins
        std::vector<ScoreScale> scales;
        std::vector<std::vector<double> > w = example.w;
        std::vector<Scores> scaleSets = example.trainBins();
        for (std::size_t b = 0; b < Example::numBins; ++b) {
            bins[b].calcScores(w[b]);
            scales.push_back(scaleSets[b].calcScoreScale(w[b], trainDecoys, fdr));
        }
        Scores merged(false);
        merged.merge(bins, fdr, false, w, trainDecoys, scales);
        EXPECT_EQ(merged.calcQvals(testFdr), count);
        // the raw scores, with none or without rescaling (--static)
        const int rawCount = directCount(example.handMerged(false), testFdr);
        EXPECT_EQ(rawCount, Scores::countMergedTargets(
                                example.testBins(), trainBins, example.w,
                                MergeNormalization::NONE, fdr, false, testFdr));
        EXPECT_EQ(rawCount, Scores::countMergedTargets(
                                example.testBins(), trainBins, example.w,
                                trainDecoys, fdr, true, testFdr));
    }
    // no bin reaches q < 0.25 on its own, the merged list does
    {
        std::vector<Scores> bins = example.testBins();
        int binsCount = 0;
        for (std::size_t b = 0; b < Example::numBins; ++b) {
            std::vector<double> w = example.w[b];
            binsCount += bins[b].calcScoresAndQvals(w, 0.25);
        }
        EXPECT_EQ(0, binsCount);
    }
    // a training set with a single decoy gives no standard deviation
    {
        std::vector<Scores> trainBins = example.trainBins();
        Scores oneDecoy(false);
        for (const ScoreHolder& sh : trainBins[1]) {
            if (sh.isTarget() || sh.pPSM->features[0] == 0.0) {
                oneDecoy.addScoreHolder(sh);
            }
        }
        oneDecoy.recalculateSizes();
        trainBins[1] = oneDecoy;
        EXPECT_EQ(-1, Scores::countMergedTargets(example.testBins(), trainBins,
                                                 example.w, trainDecoys, fdr,
                                                 false, 0.25));
    }
    // one bin: rescaling does not change the order, and the pooled count is
    // the bin's own count, with every scheme
    for (MergeNormalization scheme :
         {MergeNormalization::TEST, MergeNormalization::TRAIN, trainDecoys,
          MergeNormalization::NONE}) {
        SCOPED_TRACE(static_cast<int>(scheme));
        std::vector<Scores> bins = example.testBins(1);
        std::vector<std::vector<double> > w(1, example.w[0]);
        const int count = Scores::countMergedTargets(
            bins, example.trainBins(1), w, scheme, fdr, false, 0.3);
        EXPECT_EQ(4, count);
        EXPECT_EQ(bins[0].calcScoresAndQvals(w[0], 0.3), count);
    }
    FeatureNames::resetNumFeatures();
}

// A merged list whose targets all precede at least 200 decoys separates too
// well: the p values of all targets are below the smallest lambda of the pi0
// estimate, and the final list's postMergeStep ends the run. Two bins of 150
// targets above numDecoys decoys each: with 150 decoys per bin
// countMergedTargets returns -1 with every scheme, and merge, which gives the
// final list, throws; with 50 decoys per bin (101 p value steps) the list is
// usable and all 300 targets have q < 0.01.
TEST_F(ScoresTest, CountMergedTargetsOfATooGoodSeparationIsUnusable)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(1);
    const std::size_t numBins = 2, numTargets = 150;
    const double fdr = 0.01;
    for (std::size_t numDecoys : {std::size_t(150), std::size_t(50)}) {
        SCOPED_TRACE(numDecoys);
        const bool tooGood = numDecoys == 150u;
        const std::size_t binSize = numTargets + numDecoys;
        std::vector<double> features(numBins * binSize);
        std::vector<PSMDescription> psms(features.size());
        std::vector<Scores> bins(numBins, Scores(false));
        for (std::size_t bin = 0; bin < numBins; ++bin) {
            for (std::size_t j = 0; j < binSize; ++j) {
                const std::size_t i = bin * binSize + j;
                const bool isTarget = j < numTargets;
                // targets in [10, 11.5), decoys in (-1.5, 0], both spread
                features[i] = isTarget ? 10.0 + 0.01 * static_cast<double>(j)
                                       : -0.01 * static_cast<double>(j);
                psms[i].features = &features[i];
                psms[i].scan = static_cast<unsigned int>(i);
                bins[bin].addScoreHolder(ScoreHolder(
                    0.0, isTarget ? LabelType::TARGET : LabelType::DECOY,
                    &psms[i]));
            }
            bins[bin].recalculateSizes();
        }
        const std::vector<std::vector<double> > w(numBins, {1.0, 0.0});
        for (MergeNormalization scheme :
             {MergeNormalization::TEST, MergeNormalization::TRAIN,
              MergeNormalization::TRAIN_DECOYS, MergeNormalization::NONE}) {
            SCOPED_TRACE(static_cast<int>(scheme));
            // the bins serve as their own training sets
            EXPECT_EQ(tooGood ? -1 : 300,
                      Scores::countMergedTargets(bins, bins, w, scheme, fdr,
                                                 false, fdr));
        }
        EXPECT_EQ(tooGood ? -1 : 300,
                  Scores::countMergedTargets(bins, bins, w,
                                             MergeNormalization::NONE, fdr,
                                             true, fdr));
        std::vector<Scores> scored(bins);
        std::vector<std::vector<double> > wMerge(w);
        for (std::size_t bin = 0; bin < numBins; ++bin) {
            scored[bin].calcScores(wMerge[bin]);
        }
        Scores merged(false);
        if (tooGood) {
            EXPECT_THROW(merged.merge(scored, fdr, false, wMerge,
                                      MergeNormalization::NONE),
                         MyException);
        } else {
            merged.merge(scored, fdr, false, wMerge, MergeNormalization::NONE);
            EXPECT_EQ(300, merged.calcQvals(fdr));
        }
    }
    FeatureNames::resetNumFeatures();
}

// SanityCheck on the hand-made example. The default direction, lower x,
// separates nothing in either check. No bin reaches q < 0.25 with the
// trained weights on its own, so the bins check resets to the default
// direction; the merged list holds 12 targets with q < 0.25, so the pooled
// check keeps the trained weights and only then scores the bins with them.
// With one bin both checks give the same counts and decision.
TEST_F(ScoresTest, PooledDirectionCheckOnHandMadeBins)
{
    typedef PooledCountExample Example;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(1);
    Example example;
    const std::vector<double> defaultWeights = {-1.0, 0.0};
    for (std::size_t numBins : {Example::numBins, std::size_t(1)}) {
        const double testFdr = numBins == 1u ? 0.3 : 0.25;
        const std::vector<std::vector<double> > trained(
            example.w.begin(), example.w.begin() + numBins);
        std::vector<Scores> binsAfter[2];
        for (bool pooled : {false, true}) {
            SCOPED_TRACE(::testing::Message() << numBins << " bins, pooled "
                                              << pooled);
            SanityCheck check;
            check.addDefaultWeights(defaultWeights);
            check.setDirectionCheck(pooled
                                        ? SanityCheck::DirectionCheck::POOLED
                                        : SanityCheck::DirectionCheck::BINS);
            check.setMergeScheme(MergeNormalization::TRAIN_DECOYS, 0.01, false);
            std::vector<Scores> bins = example.testBins(numBins);
            std::vector<Scores> trainBins = example.trainBins(numBins);
            std::vector<std::vector<double> > w(numBins,
                                                std::vector<double>(2));
            EXPECT_EQ(0, check.getInitDirection(bins, trainBins, NULL, w,
                                                testFdr, testFdr));
            EXPECT_EQ(std::vector<std::vector<double> >(numBins,
                                                        defaultWeights),
                      w);
            w = trained;
            const bool keep = numBins == 1u || pooled;
            EXPECT_EQ(keep, check.validateDirection(bins, trainBins, w));
            if (keep) {
                EXPECT_EQ(trained, w);
                // the bins are scored with their folds' weights
                for (std::size_t b = 0; b < numBins; ++b) {
                    for (const ScoreHolder& sh : bins[b]) {
                        EXPECT_DOUBLE_EQ(
                            w[b][0] * sh.pPSM->features[0] + w[b][1],
                            sh.score);
                    }
                }
            } else {
                EXPECT_EQ(std::vector<std::vector<double> >(numBins,
                                                            defaultWeights),
                          w);
            }
            binsAfter[pooled] = bins;
        }
        if (numBins == 1u) {
            ASSERT_EQ(binsAfter[0][0].size(), binsAfter[1][0].size());
            for (auto it = binsAfter[0][0].begin(),
                      jt = binsAfter[1][0].begin();
                 it != binsAfter[0][0].end(); ++it, ++jt) {
                EXPECT_EQ(it->pPSM, jt->pPSM);
                EXPECT_EQ(it->score, jt->score);
                EXPECT_EQ(it->q, jt->q);
            }
        }
    }
    FeatureNames::resetNumFeatures();
}

// The second chance on the hand-made example, whose default direction,
// lower x, separates nothing. At q < 0.25 no bin reaches the test FDR on
// its own, so the per-bin count would reset; the merged list rescaled by
// the training decoys holds 12 targets with q < 0.25, the initial
// direction none, so the trained weights are kept, the bins stay scored
// with them and the merge is to use train-decoys. At q < 0.05, which the
// first 12 targets (q = 1/12) miss, the merged list holds none and the
// check resets as bins does. With --override the decision is that of bins,
// and with one bin there is no second chance. The check leaves the random
// number generator as it found it.
TEST_F(ScoresTest, SecondChanceOnHandMadeBins)
{
    typedef PooledCountExample Example;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(1);
    Example example;
    const std::vector<double> defaultWeights = {-1.0, 0.0};
    struct Case {
        std::size_t numBins;
        double testFdr;
        bool overrule;
        bool keep;
    };
    const Case cases[] = {{Example::numBins, 0.25, false, true},
                          {Example::numBins, 0.05, false, false},
                          {Example::numBins, 0.25, true, false},
                          {1u, 0.25, false, false}};
    for (const Case& c : cases) {
        SCOPED_TRACE(::testing::Message() << c.numBins << " bins, testFdr "
                                          << c.testFdr << ", override "
                                          << c.overrule);
        SanityCheck check;
        check.addDefaultWeights(defaultWeights);
        check.setDirectionCheck(SanityCheck::DirectionCheck::SECOND_CHANCE);
        check.setOverrule(c.overrule);
        // the scheme of the pooled check; the second chance counts with
        // train-decoys whatever it is
        check.setMergeScheme(MergeNormalization::TEST, 0.01, false);
        std::vector<Scores> bins = example.testBins(c.numBins);
        std::vector<Scores> trainBins = example.trainBins(c.numBins);
        std::vector<std::vector<double> > w(c.numBins,
                                            std::vector<double>(2));
        EXPECT_EQ(0, check.getInitDirection(bins, trainBins, NULL, w,
                                            c.testFdr, c.testFdr));
        const std::vector<std::vector<double> > trained(
            example.w.begin(), example.w.begin() + c.numBins);
        w = trained;
        PseudoRandom::setSeed(7u);
        const uint64_t randomState = PseudoRandom::getSeed();
        EXPECT_EQ(c.keep, check.validateDirection(bins, trainBins, w));
        EXPECT_EQ(randomState, PseudoRandom::getSeed());
        EXPECT_EQ(c.keep, check.secondChanceKept());
        EXPECT_EQ(MergeNormalization::TRAIN_DECOYS,
                  check.secondChanceMergeScheme());
        if (c.keep || c.overrule) {
            EXPECT_EQ(trained, w);
        } else {
            EXPECT_EQ(std::vector<std::vector<double> >(c.numBins,
                                                        defaultWeights),
                      w);
        }
        if (c.keep) {
            // the per-bin count scored the bins with their folds' weights
            for (std::size_t b = 0; b < c.numBins; ++b) {
                for (const ScoreHolder& sh : bins[b]) {
                    EXPECT_DOUBLE_EQ(w[b][0] * sh.pPSM->features[0] + w[b][1],
                                     sh.score);
                }
            }
        }
    }
    FeatureNames::resetNumFeatures();
}

// Two test bins of 6 targets and 3 decoys, with two features: f0 puts all
// targets of a bin above its decoys, but the second bin's f0 is shifted by
// +10 (targets 20..15, decoys 11, 10, 9), which its training set does not
// show; f1 is not shifted and puts 3 targets of each bin above its first
// decoy. The training sets choose f0 as the default direction. At
// --testFDR 0.2 the trained vector f1 fails the per-bin count (q = 1/3 at
// best), so the second chance runs:
//   trained f1, merged with train-decoys:          6 targets with q < 0.2
//   default f0, merged with train-decoys:          6 (the shifted bin's
//                                                  decoys precede the
//                                                  other bin's targets)
//   default f0, merged with test, as reported:    12 (each bin anchored at
//                                                  its own scores)
// The reset, which merges with the run's test scheme, would report 12
// targets, so the second chance must not keep the trained vector, as a
// comparison with the initial direction counted with train-decoys (6) did.
// With --init-weights f1 the initial direction is f1, but a reset gives the
// default vectors f0, which the second chance therefore counts.
class ResetCountExample {
  public:
    static constexpr std::size_t numBins = 2;
    ResetCountExample() {
        const double tF0[6] = {10, 9, 8, 7, 6, 5}, dF0[3] = {1, 0, -1};
        const double tF1[6] = {5, 4, 3, 2, 1, -2}, dF1[3] = {2.5, 0, -1};
        // training sets: f0 separates all 4 targets, f1 only 2
        const double trF0[7] = {4, 5, 6, 7, -1, 0, 1};
        const double trF1[7] = {0.5, 1.5, -0.5, 2, -1, 0, 1};
        for (std::size_t bin = 0; bin < numBins; ++bin) {
            const double shift = bin == 1u ? 10.0 : 0.0;
            for (std::size_t i = 0; i < 6; ++i) {
                add(test[bin], tF0[i] + shift, tF1[i], true);
            }
            for (std::size_t i = 0; i < 3; ++i) {
                add(test[bin], dF0[i] + shift, dF1[i], false);
            }
            for (std::size_t i = 0; i < 7; ++i) {
                add(train[bin], trF0[i], trF1[i], i < 4u);
            }
        }
    }
    std::vector<Scores> testBins() const { return sets(test); }
    std::vector<Scores> trainBins() const { return sets(train); }
    const std::vector<std::vector<double> > f0 = {{1, 0, 0}, {1, 0, 0}};
    const std::vector<std::vector<double> > f1 = {{0, 1, 0}, {0, 1, 0}};

  private:
    struct Psm {
        double features[2];
        PSMDescription psm;
        bool isTarget;
    };
    void add(std::vector<std::unique_ptr<Psm> >& set, double x0, double x1,
             bool isTarget) {
        set.emplace_back(new Psm());
        Psm& p = *set.back();
        p.features[0] = x0;
        p.features[1] = x1;
        p.psm.features = p.features;
        p.psm.scan = static_cast<unsigned int>(++numPsms);
        p.isTarget = isTarget;
    }
    std::vector<Scores> sets(
        const std::vector<std::unique_ptr<Psm> > (&psms)[numBins]) const {
        std::vector<Scores> result(numBins, Scores(false));
        for (std::size_t bin = 0; bin < numBins; ++bin) {
            for (const std::unique_ptr<Psm>& p : psms[bin]) {
                result[bin].addScoreHolder(ScoreHolder(
                    0.0, p->isTarget ? LabelType::TARGET : LabelType::DECOY,
                    &p->psm));
            }
            result[bin].recalculateSizes();
        }
        return result;
    }
    std::vector<std::unique_ptr<Psm> > test[numBins], train[numBins];
    std::size_t numPsms = 0;
};

TEST_F(ScoresTest, SecondChanceComparesWithTheReportedReset)
{
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(2);
    ResetCountExample example;
    const double testFdr = 0.2, fdr = 0.01;
    typedef MergeNormalization Scheme;
    // the three counts of the comment above
    EXPECT_EQ(6, Scores::countMergedTargets(
                     example.testBins(), example.trainBins(), example.f1,
                     Scheme::TRAIN_DECOYS, fdr, false, testFdr));
    EXPECT_EQ(6, Scores::countMergedTargets(
                     example.testBins(), example.trainBins(), example.f0,
                     Scheme::TRAIN_DECOYS, fdr, false, testFdr));
    EXPECT_EQ(12, Scores::countMergedTargets(
                      example.testBins(), example.trainBins(), example.f0,
                      Scheme::TEST, fdr, false, testFdr));
    // --init-weights is read through the normalizer, the identity here
    std::unique_ptr<Normalizer> norm = Normalizer::create(Normalizer::NONORM);
    norm->setNumFeatures(2);
    const std::string weightFile =
        ::testing::TempDir() + "SecondChanceComparesWithTheReportedReset.txt";
    {
        std::ofstream out(weightFile.c_str());
        out << "f0\tf1\tm0\n0\t1\t0\n0\t1\t0\n";
    }
    for (bool initWeights : {false, true}) {
        // with train-decoys as the run's scheme, as with more than three
        // folds, the reset reports 6 targets and the trained vector is kept
        for (Scheme scheme : {Scheme::TEST, Scheme::TRAIN_DECOYS}) {
            SCOPED_TRACE(::testing::Message() << "--init-weights "
                                              << initWeights << ", scheme "
                                              << static_cast<int>(scheme));
            SanityCheck check;
            check.setDirectionCheck(SanityCheck::DirectionCheck::SECOND_CHANCE);
            check.setMergeScheme(scheme, fdr, false);
            if (initWeights) {
                check.setInitWeightFN(weightFile);
            }
            std::vector<Scores> bins = example.testBins();
            std::vector<Scores> trainBins = example.trainBins();
            std::vector<std::vector<double> > w(
                ResetCountExample::numBins, std::vector<double>(3));
            // f0 per bin: 6 targets above 3 decoys, q = 1/6; f1: none
            EXPECT_EQ(initWeights ? 0 : 12,
                      check.getInitDirection(bins, trainBins, norm.get(), w,
                                             testFdr, testFdr));
            EXPECT_EQ(initWeights ? example.f1 : example.f0, w);
            w = example.f1;
            PseudoRandom::setSeed(11u);
            const bool keep = scheme == Scheme::TRAIN_DECOYS;
            EXPECT_EQ(keep, check.validateDirection(bins, trainBins, w));
            EXPECT_EQ(11u, PseudoRandom::getSeed());
            EXPECT_EQ(keep, check.secondChanceKept());
            // kept, or reset to the default vectors, not to --init-weights
            EXPECT_EQ(keep ? example.f1 : example.f0, w);
        }
    }
    std::remove(weightFile.c_str());
    FeatureNames::resetNumFeatures();
}

// countTargetsBelowFdr must report the same number of targets below the FDR
// threshold as calcScoresAndQvals, also with many tied scores, with and
// without the mix-max correction.
TEST(ScoresCountTest, CountTargetsBelowFdrMatchesCalcScoresAndQvals)
{
    const size_t numFeatures = 2;
    FeatureNames::resetNumFeatures();
    FeatureNames::setNumFeatures(numFeatures);
    const int numPsms = 2000;
    std::vector<double> features(numPsms * numFeatures);
    std::vector<PSMDescription*> psms;
    for (bool usePi0 : {false, true}) {
        Scores scores(usePi0);
        unsigned int seed = 12345u;
        for (int i = 0; i < numPsms; ++i) {
            seed = seed * 1103515245u + 12345u;
            bool isTarget = (seed >> 16) % 2 == 0;
            double* row = &features[i * numFeatures];
            // coarse values, so that many PSMs get exactly the same score
            row[0] = static_cast<double>((seed >> 8) % 7) + (isTarget ? 3 : 0);
            row[1] = static_cast<double>((seed >> 4) % 3);
            PSMDescription* psm = new PSMDescription();
            psm->features = row;
            psm->scan = i / 2;
            psms.push_back(psm);
            scores.addScoreHolder(ScoreHolder(0.0,
                isTarget ? LabelType::TARGET : LabelType::DECOY, psm));
        }
        scores.recalculateSizes();
        std::vector<double> w = {1.0, -0.5, 0.25};
        scores.calcScores(w);
        scores.postMergeStep();  // sets pi0 when usePi0 is true
        if (usePi0) {
            ASSERT_LT(scores.getPi0(), 1.0);
        }
        Scores::FeatureRows rows;
        scores.getFeatureRows(rows);
        std::vector<std::pair<double, bool> > combined;
        // countTargetsBelowFdr needs the q values of the decoys as well;
        // calcScoresAndQvals sets this only after the first call below
        PosteriorEstimator::setNegative(true);
        for (double fdr : {0.01, 0.05, 0.2}) {
            for (bool skipDecoysPlusOne : {false, true}) {
                int count = scores.countTargetsBelowFdr(
                    rows, w, fdr, skipDecoysPlusOne, combined);
                EXPECT_EQ(scores.calcScoresAndQvals(w, fdr, skipDecoysPlusOne),
                          count);
                if (fdr > 0.1) {
                    EXPECT_GT(count, 0);  // the test needs accepted targets
                }
            }
        }
    }
    for (PSMDescription* psm : psms) {
        delete psm;
    }
    FeatureNames::resetNumFeatures();
}

// Both PEP pathways of calcPep, from target q values (default) and from the
// target-decoy labels (tdc), with I-spline and PAVA regression, give PEPs in
// [0, 1] that do not decrease down the list, and whose mean over the targets
// above a q value threshold is close to that threshold.
TEST(ScoresPepTest, BothPathwaysGiveMonotoneCalibratedPeps)
{
    const int numPsms = 6000;
    std::vector<PSMDescription> psms(numPsms);
    unsigned int seed = 777u;
    auto uniform = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return ((seed >> 8) + 0.5) / 16777216.0;
    };
    auto normal = [&uniform]() {
        return std::sqrt(-2.0 * std::log(uniform())) *
               std::cos(6.283185307179586 * uniform());
    };
    for (bool tdc : {false, true}) {
        for (bool pava : {false, true}) {
            Scores scores(false);
            seed = 777u;
            for (int i = 0; i < numPsms; ++i) {
                psms[i].scan = static_cast<unsigned int>(i);
                const bool isTarget = uniform() < 0.5;
                const bool correct = isTarget && uniform() < 0.4;
                scores.addScoreHolder(ScoreHolder(normal() + (correct ? 3.0 : 0.0),
                    isTarget ? LabelType::TARGET : LabelType::DECOY, &psms[i]));
            }
            scores.recalculateSizes();
            // calcQvals and calcPep expect the PSMs sorted by score
            std::sort(scores.begin(), scores.end(), std::greater<ScoreHolder>());
            scores.calcQvals(0.01);
            scores.calcPep(false, true, pava, tdc);
            double lastTargetPep = 0.0, sum = 0.0;
            int n = 0;
            for (const ScoreHolder& sh : scores) {
                ASSERT_GE(sh.pep, 0.0);
                ASSERT_LE(sh.pep, 1.0);
                if (!sh.isTarget()) continue;
                EXPECT_GE(sh.pep, lastTargetPep - 1e-12);
                lastTargetPep = sh.pep;
                if (sh.q <= 0.05) {
                    sum += sh.pep;
                    ++n;
                }
            }
            ASSERT_GT(n, 100);
            EXPECT_GT(sum / n, 0.05 / 2.0) << "tdc=" << tdc << " pava=" << pava;
            EXPECT_LT(sum / n, 0.05 * 2.0) << "tdc=" << tdc << " pava=" << pava;
        }
    }
}
