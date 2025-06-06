#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "CompositionSorter.h"
#include "Globals.h"
#include "Scores.h"

// helper for sorting ScoreHolder references by scan and score, to help TDC
bool compareByScanThenScore(const ScoreHolder* lhs, const ScoreHolder* rhs) {
  // sort by file, scan, mass in ascending order
  if (lhs->pPSM->specFileNr != rhs->pPSM->specFileNr) {
    return lhs->pPSM->specFileNr < rhs->pPSM->specFileNr;
  }
  if (lhs->pPSM->scan != rhs->pPSM->scan) {
    return lhs->pPSM->scan < rhs->pPSM->scan;
  }
  if (lhs->pPSM->expMass != rhs->pPSM->expMass) {
    return lhs->pPSM->expMass < rhs->pPSM->expMass;
  }
  // Then sort by score in descending order
  return lhs->score > rhs->score;
}

void targetDecoyCompetition(std::vector<ScoreHolder*>& scoreHolders) {
  std::sort(scoreHolders.begin(), scoreHolders.end(), compareByScanThenScore);
  if (VERB > 1) {
    std::cerr << "Before TDC there are " << scoreHolders.size() << " PSMs."
              << std::endl;
  }
  // First, sort the vector using the custom comparator

  // Use a map to track scan numbers and keep only the best score for each scan
  std::map<int, ScoreHolder*> bestScoreByScan;
  for (ScoreHolder* holder : scoreHolders) {
    int scan = holder->pPSM->scan;
    if (bestScoreByScan.find(scan) == bestScoreByScan.end() ||
        bestScoreByScan[scan]->score < holder->score) {
      bestScoreByScan[scan] = holder;
    }
  }

  // Clear the original scoreHolders and repopulate it with the best scores from
  // each scan
  scoreHolders.clear();
  for (const auto& pair : bestScoreByScan) {
    scoreHolders.push_back(pair.second);
  }
  if (VERB > 1) {
    std::cerr << "After TDC there are " << scoreHolders.size() << " PSMs."
              << std::endl;
  }
}

std::string CompositionSorter::generateCompositionSignature(
    const std::string& peptide) {
  std::map<std::string, int> counts;

  for (size_t i = 0; i < peptide.size(); ++i) {
    if (peptide[i] == '[') {
      size_t j = i + 1;
      while (j < peptide.size() && peptide[j] != ']') {
        ++j;
      }
      // Extract the modification
      counts[peptide.substr(i, j - i + 1)]++;
      // Skip the modification
      i = j;
    } else {
      counts[std::string(1, peptide[i])]++;
    }
  }

  std::string signature;
  for (const auto& aminoAcidCountPair : counts) {
    const auto& aminoAcid = aminoAcidCountPair.first;
    const auto& count = aminoAcidCountPair.second;
    signature += aminoAcid;
    signature += std::to_string(count);
  }
  return signature;
}

int CompositionSorter::addPSMs(const Scores& scores, bool useTDC) {
  if (useTDC) {
    std::vector<ScoreHolder*> scoreHolders;
    for (const auto& scr : scores) {
      scoreHolders.push_back(const_cast<ScoreHolder*>(&scr));
    }
    targetDecoyCompetition(scoreHolders);
    for (auto& pScr : scoreHolders) {
      std::string peptide = pScr->getPSM()->getPeptideSequence();
      std::string signature = generateCompositionSignature(peptide);
      compositionToPeptidesToScore_[signature][peptide].push_back(pScr);
    }
  } else {
    for (const auto& scr : scores) {
      std::string peptide = scr.getPSM()->getPeptideSequence();
      std::string signature = generateCompositionSignature(peptide);
      compositionToPeptidesToScore_[signature][peptide].push_back(
          const_cast<ScoreHolder*>(&scr));
    }
  }
  return 0;
}

int CompositionSorter::sortScorePerPeptide() {
  // Comparator for sorting ScoreHolder references
  auto compareScoreHolder = [](const ScoreHolder* lhs, const ScoreHolder* rhs) {
    return lhs->score > rhs->score;  // For descending order
  };

  // Iterating and sorting
  for (auto& compositionToPeptides : compositionToPeptidesToScore_) {
    auto& composition = compositionToPeptides.first;
    auto& peptideMap = compositionToPeptides.second;

    for (auto& peptideToScores : peptideMap) {
      auto& peptide = peptideToScores.first;
      auto& scores = peptideToScores.second;

      std::sort(scores.begin(), scores.end(), compareScoreHolder);
    }
  }
  return 0;
}

void incVector(std::vector<size_t>& sizes, size_t s) {
  // Ensure the vector is large enough to hold the element at index 's'
  if (s >= sizes.size()) {
    // Resize the vector to accommodate the new index, initializing new elements
    // to 0
    sizes.resize(s + 1, 0);
  }
  // Increment the element at index 's'
  sizes[s]++;
}

int CompositionSorter::inCompositionCompetition(Scores& bestScoreHolders,
                                                unsigned int decoysPerTarget) {
  // Sort the PSMs for each peptide in descending order of score
  sortScorePerPeptide();

  if (VERB > 1) {
    std::cerr << "Composition Matching starting with "
              << compositionToPeptidesToScore_.size() << " compositions."
              << std::endl;
  }

  std::vector<size_t> compSizeStat;
  std::vector<size_t> compTargetSizeStat;
  std::vector<size_t> compNumTargetSizeStat;
  for (auto& compositionToPeptides : compositionToPeptidesToScore_) {
    auto& composition = compositionToPeptides.first;
    auto& peptideMap = compositionToPeptides.second;
    std::vector<std::vector<ScoreHolder*>> compositionGroups;
    std::vector<ScoreHolder*> targets;
    std::vector<ScoreHolder*> decoys;
    // register length statistics for printouts
    size_t numberOfPeptidesInComposition = peptideMap.size();
    incVector(compSizeStat, numberOfPeptidesInComposition);
    // Add the target peptides
    size_t numTargetPeptides = 0;
    for (auto& peptideToScores : peptideMap) {
      auto& peptide = peptideToScores.first;
      auto& scoreHolders = peptideToScores.second;
      if (scoreHolders.empty()) {
        continue;
      }
      ScoreHolder* firstScoreHolder =
          scoreHolders.front();  // This is the best scoring PSM for the peptide
      if (firstScoreHolder->isTarget()) {
        targets.push_back(firstScoreHolder);
        numTargetPeptides++;
        incVector(compNumTargetSizeStat, numberOfPeptidesInComposition);
      } else {
        decoys.push_back(firstScoreHolder);
      }
    }
    if (numTargetPeptides > 0)
      incVector(compTargetSizeStat, numberOfPeptidesInComposition);
    // Format tuples of targets and decoys. Now we are not checking if there are
    // enough decoys for each target.
    // TODO: Check if there are enough decoys for each target
    std::reverse(decoys.begin(), decoys.end());
    for (ScoreHolder* target : targets) {
      std::vector<ScoreHolder*> group;
      group.push_back(target);  // Add the target ScoreHolder to the group

      // Add decoysPerTarget number of decoys to the group
      for (unsigned int i = 0; i < decoysPerTarget && !decoys.empty(); ++i) {
        group.push_back(decoys.back());  // Add the last decoy (the one first
                                         // added to the vector)
        decoys.pop_back();  // Remove the added decoy from the decoys list
      }
      compositionGroups.push_back(
          group);  // Add the newly formed group to compositionGroups
    }
    // Take care of the remaining decoys
    for (ScoreHolder* decoy : decoys) {
      std::vector<ScoreHolder*> group;
      group.push_back(decoy);  // Add the target ScoreHolder to the group

      // Add decoysPerTarget number of decoys to the group
      for (unsigned int i = 1; i < decoysPerTarget && !decoys.empty(); ++i) {
        group.push_back(decoys.back());  // Add the last decoy (the one first
                                         // added to the vector)
        decoys.pop_back();  // Remove the added decoy from the decoys list
      }
      compositionGroups.push_back(
          group);  // Add the newly formed group to compositionGroups
    }
    // Sort each group in compositionGroups based on the ScoreHolder's score
    // and add it to bestScoreHolders
    for (auto& group : compositionGroups) {
      if (group.empty())
        continue;
      std::sort(group.begin(), group.end(),
                [](const ScoreHolder* a, const ScoreHolder* b) -> bool {
                  return a->score >
                         b->score;  // Sort in descending order of score
                });
      bestScoreHolders.addScoreHolder(
          *(group.front()));  // Add the highest-scoring ScoreHolder
    }
  }
  if (VERB > 1) {
    std::cerr
        << "Composition Group Statistics: Sizes of composition groups, sizes "
           "with at least 1 target peptide, and avg. num targets per group."
        << std::endl;
    std::cerr << "Size\tTotal\tIsTarget\tAvgTarget" << std::endl;
    for (size_t ix = compSizeStat.size(); --ix;) {
      std::cerr << ix << '\t' << compSizeStat[ix] << '\t';
      if (ix < compTargetSizeStat.size()) {
        std::cerr << compTargetSizeStat[ix] << '\t'
                  << compTargetSizeStat[ix] / ix << std::endl;
      } else {
        std::cerr << 0 << '\t' << 0 << std::endl;
      }
    }
  }

  int numTargets = 0;
  for (const ScoreHolder& pSH : bestScoreHolders) {
    if (pSH.isTarget())
      numTargets++;
  }
  if (VERB > 1) {
    std::cerr << "Composition Matching ends with " << bestScoreHolders.size()
              << " peptides, whereof " << numTargets << " is target peptides."
              << std::endl;
  }
  return 0;
}

int CompositionSorter::psmAndPeptide(const Scores& scores,
                                     Scores& winnerPeptides,
                                     unsigned int decoysPerTarget) {
  //    psmLevelCompetition(); // Should be handled by the reader?

  // Populate the structure with the winning PSMs
  addPSMs(scores);

  // Split out tuples of peptides of identical composition, and select the most
  // high scoring peptide in each tuple
  inCompositionCompetition(winnerPeptides, decoysPerTarget);
  winnerPeptides.recalculateSizes();
  return 0;
}

void CompositionSorter::psmsOnly(const Scores& scores, Scores& winnerPeptides) {
  // Use an unordered_map (hash map) to store the best ScoreHolder for each
  // peptide sequence
  std::unordered_map<std::string, ScoreHolder*> bestScoreHolders;

  // Iterate over each ScoreHolder
  for (const ScoreHolder& sh : scores) {
    const std::string& peptide = sh.getPSM()->getPeptideSequence();
    // Update the map if the current ScoreHolder has a higher score for the
    // peptide
    if (bestScoreHolders.find(peptide) == bestScoreHolders.end() ||
        sh.score > bestScoreHolders[peptide]->score) {
      bestScoreHolders[peptide] = const_cast<ScoreHolder*>(&sh);
    }
  }

  // Collect the best ScoreHolders into a result vector
  for (const auto& entry : bestScoreHolders) {
    winnerPeptides.addScoreHolder(*(entry.second));
  }
  winnerPeptides.recalculateSizes();

  if (VERB > 1) {
    std::cerr << "Selected the best scoring PSM for each of the "
              << winnerPeptides.size() << " peptides from a dataset of "
              << scores.size() << " PSMs." << std::endl;
  }
}

void CompositionSorter::retainRepresentatives(const Scores& psms,
                                              Scores& winnerPeptides,
                                              double selectionFDR,
                                              unsigned int decoysPerTarget,
                                              bool useCompositionMatch) {
  if (useCompositionMatch) {
    if (VERB > 1) {
      std::cerr << "Starting reset: psmAndPeptide" << std::endl;
    }
    CompositionSorter sorter;
    sorter.psmAndPeptide(psms, winnerPeptides, decoysPerTarget);
  } else {
    if (VERB > 1) {
      std::cerr << "Starting reset: psmsOnly" << std::endl;
    }
    psmsOnly(psms, winnerPeptides);
  }
}