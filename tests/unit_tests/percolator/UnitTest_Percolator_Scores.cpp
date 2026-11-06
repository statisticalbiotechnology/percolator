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
#include <cstdarg>
#include <map>
#include <random>
#include <utility>
#include "SetHandler.h"
#include "DataSet.h"
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
